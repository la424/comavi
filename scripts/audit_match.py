"""Boundary-aware literal matching for the manuscript audits.

WHY THIS MODULE EXISTS
======================
A plain ``needle in text`` test cannot see a value that prose has made MORE
precise than the data supports. ``"0.903" in "0.9039"`` is True, so the
manuscript could report an AUC of 0.9039 -- a number no analysis produced --
and an audit gating on the literal ``0.903`` would still pass.

This is not hypothetical. Mutation-testing every digit-ending gate in
``audit_evidence_claims.py`` and ``audit_tier_energy_claims.py`` (append one
digit to the value where the prose states it, then re-run the audit) found
**45 of 45 gates non-binding**, including both paired discrimination AUCs, the
exactly enumerated permutation p, and every screen confidence bound. The same
defect had already been found and fixed once at a single call site; this module
exists so the fix is structural rather than per-literal.

Both audits route their containment tests through :func:`contains_literal`.
:func:`weak_needle` additionally refuses needles that cannot constrain anything
even with correct boundaries -- a bare ``"11"`` matches inside ``"11/14"``,
``"0.115"`` and a dozen unrelated places, so it must carry context.

Run ``python scripts/audit_match.py`` to execute this module's self-tests.
"""
from __future__ import annotations


def _extends_left(text: str, i: int) -> bool:
    """True if the match at ``i`` is the tail of a longer number."""
    if i == 0:
        return False
    if text[i - 1].isdigit():
        return True
    # ".24" inside "0.24": needle is the fractional tail of a decimal.
    return text[i - 1] == "." and i > 1 and text[i - 2].isdigit()


def _extends_right(text: str, j: int) -> bool:
    """True if the match ending at ``j`` is the head of a longer number."""
    if j >= len(text):
        return False
    if text[j].isdigit():
        return True
    # "0" inside "0.24": needle is the integer part of a decimal.
    return text[j] == "." and j + 1 < len(text) and text[j + 1].isdigit()


def contains_literal(text: str, needle: str) -> bool:
    """``needle in text``, but a numeric needle must match a COMPLETE number.

    Boundary rules are applied only at the ends where the needle is numeric,
    so phrase needles ("the 25 gradeable E1 axes") behave as before while
    numeric needles ("0.903", "17/30", "p = 0.0095") cannot be satisfied by a
    longer, more precise number.
    """
    # A needle is numeric at an edge if it starts/ends with a digit, or with a
    # decimal point adjacent to one (".24" is the tail of "0.24", so its left
    # edge needs the same guard as a bare digit would).
    guard_left = needle[:1].isdigit() or (needle[:1] == "." and needle[1:2].isdigit())
    guard_right = needle[-1:].isdigit() or (needle[-1:] == "." and needle[-2:-1].isdigit())
    start = 0
    while True:
        i = text.find(needle, start)
        if i < 0:
            return False
        j = i + len(needle)
        if not (guard_left and _extends_left(text, i)) and \
           not (guard_right and _extends_right(text, j)):
            return True
        start = i + 1


def check_literal(text: str, needle: str) -> str | None:
    """Full gate for one literal: return an error message, or None if clean.

    Two failure modes, both observed in this project:

    1. **Absent** -- no complete-number occurrence of the value.
    2. **Partially stale** -- the value appears correctly somewhere AND appears
       elsewhere as the truncated head of a longer number. A presence-only test
       is satisfied by the good occurrence and cannot see the bad one. This is
       the shape of the drift that left a stale 4,000-draw bound in the prose
       after the estimator moved to 200,000 draws: one mention was updated and
       another was not.

    Mode 2 is only checked for needles that are numeric at their right edge,
    since only those can be extended into a different value.
    """
    if not contains_literal(text, needle):
        return f"absent (no complete-number match for {needle!r})"

    guard_right = needle[-1:].isdigit() or (needle[-1:] == "." and needle[-2:-1].isdigit())
    if guard_right:
        start, extended = 0, []
        while True:
            i = text.find(needle, start)
            if i < 0:
                break
            if _extends_right(text, i + len(needle)):
                extended.append(text[i:i + len(needle) + 6])
            start = i + 1
        if extended:
            return (f"{needle!r} also appears as the head of a longer number "
                    f"({', '.join(repr(e) for e in extended[:3])}) -- one mention "
                    "may have been updated while another went stale")
    return None


def weak_needle(needle: str) -> str | None:
    """Reason a needle cannot meaningfully gate a value, or None if it can.

    A short bare number carries no context: it will match somewhere in any
    manuscript regardless of what the prose claims, so a gate built on it is
    decoration. Such needles must be given surrounding words or a unit.
    """
    bare = needle.strip()
    if all(c.isdigit() or c in ".," for c in bare) and len(bare.replace(".", "")) < 3:
        return (f"bare numeric needle {needle!r} is too short to locate a claim; "
                "give it surrounding prose (e.g. 'none of the 17' not '0')")
    return None


def _self_test() -> None:
    ok = 0

    # The exact defect this module exists to prevent.
    assert "0.903" in "AUC 0.9039 on the set"                 # plain containment: fooled
    assert not contains_literal("AUC 0.9039 on the set", "0.903")
    assert contains_literal("AUC 0.903 on the set", "0.903")
    ok += 1

    # Fractions and signed / prefixed values.
    assert not contains_literal("recall 17/170", "17/17")
    assert contains_literal("recall 17/17 (perfect)", "17/17")
    assert contains_literal("CI +0.345 to", "+0.345")
    assert contains_literal("bound of 0.24]", "0.24")
    ok += 1

    # A needle preceded by a digit is the tail of a longer number.
    assert not contains_literal("p = 10.0095", "0.0095")
    assert contains_literal("p = 0.0095 (exact)", "0.0095")
    ok += 1

    # Decimal-context rules: "0" and ".24" must not match inside "0.24".
    assert not contains_literal("gain of 0.24 overall", "0")
    assert not contains_literal("gain of 0.24 overall", ".24")
    assert contains_literal("exactly 0 variants", "0")
    ok += 1

    # Phrase needles keep plain-containment behaviour at their non-numeric end.
    assert contains_literal("the 25 gradeable E1 axes are", "the 25 gradeable E1 axes")
    assert contains_literal("(K_D, Tm; n = 29) and", "K_D, Tm; n = 29)")
    ok += 1

    # Embedded-phrase needle whose numeric tail is extended.
    assert not contains_literal("agreement is 18/25 = 0.7201", "agreement is 18/25 = 0.720")
    assert contains_literal("agreement is 18/25 = 0.720.", "agreement is 18/25 = 0.720")
    ok += 1

    # check_literal additionally catches the PARTIAL-UPDATE case: value correct
    # in one mention, silently extended in another. Presence-only tests miss it.
    two_ok = "specificity 0.567 ... and 0.567 again"
    one_stale = "specificity 0.567 ... and 0.5671 again"
    assert contains_literal(one_stale, "0.567")     # presence test: fooled
    assert check_literal(two_ok, "0.567") is None
    assert check_literal(one_stale, "0.567") is not None
    assert "stale" in check_literal(one_stale, "0.567")
    assert "absent" in check_literal("no such value here", "0.567")
    # A non-numeric right edge cannot be extended, so mode 2 does not apply.
    assert check_literal("the 25 gradeable E1 axes and E1 axes", "gradeable E1 axes") is None
    ok += 1

    # weak_needle flags exactly the needles that cannot constrain a claim.
    assert weak_needle("0") is not None
    assert weak_needle("11") is not None
    assert weak_needle("0.903") is None
    assert weak_needle("17/17") is None
    assert weak_needle("none of the 17") is None
    ok += 1

    print(f"audit_match self-tests: {ok}/8 groups PASS")


if __name__ == "__main__":
    _self_test()
