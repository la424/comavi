#!/usr/bin/env python3
"""Does each committed expectation agree with the evidence written beside it?

WHY THIS EXISTS
---------------
Every expected token in COMAVI_evidence_ledger.csv carries a one-line
`evidence_basis` and a citation. Nothing checked that the token and the basis
say the same thing, and three of them did not:

  msh2_msh6 C697F  binding  destab  basis: "...protein stable, interaction intact"
  kras_craf G12D   binding  destab  basis: "WT-like RAF1 binding; impaired GTP
                                            hydrolysis. G12 is 13.6 A from RAF1"
  kras_craf G12V   binding  destab  basis: "WT-like RAF1 binding; impaired GTP
                                            hydrolysis"

All three were introduced by one v7.6 review batch that changed `neutral` to
`destab` without changing the basis. They are reversed in v7.9. This gate makes
that class of error fail loudly instead of propagating into the ground truth,
the grading, the prioritization population and the manuscript.

WHY IT IS AXIS-AWARE
--------------------
A naive "token says destab but the basis says intact" rule catches only one of
the three. The two KRAS bases contain BOTH intactness language ("WT-like RAF1
binding") and loss language ("impaired GTP hydrolysis") -- because the variant
really does break something, just not the thing the axis is about. The loss is
catalytic; the axis is binding.

So the gate matches intactness and loss assertions against the vocabulary of
the axis being committed. "WT-like RAF1 binding" is an intactness claim about
binding and conflicts with a destab binding token. "Impaired GTP hydrolysis" is
a loss claim about catalysis and is evidence for neither.

WHAT IT CANNOT DO
-----------------
It reads one curated English sentence, so it catches contradictions that are
stated, not omissions or wrong citations. A basis that says nothing about the
axis is not flagged -- silence is not contradiction. Treat a pass as "no token
contradicts its own written evidence", not as "every token is correct".

Usage
-----
    python scripts/verify_ledger_token_basis.py [--ledger PATH] [--verbose]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "reference_outputs" / "COMAVI_evidence_ledger.csv"

# Vocabulary that identifies which axis a clause is talking about. Partner and
# gene names are deliberately included: "WT-like RAF1 binding" and "abolished
# phospho-Smad2/3 interaction" are both binding clauses.
AXIS_VOCAB = {
    "binding": r"(bind\w*|interact\w*|affinit\w*|complex\w*|partner|assembl\w*|"
               r"adhesion|heterodim\w*|K_?D|dissociation|co-?IP|two-?hybrid|SPR)",
    "monomer": r"(fold\w*|stabilit\w*|stable|unfold\w*|denatur\w*|Tm\b|thermal|"
               r"subunit|monomer|native\w*|CD\b|DSC|abundance|expression|degrad\w*)",
    "fold_complex": r"(fold\w*|stabilit\w*|stable|unfold\w*|denatur\w*|Tm\b|thermal|"
                    r"assembl\w*|complex\w*|native\w*)",
}

# Stemmed throughout: an earlier version matched "retained" but not "retains",
# and read "Retains p110 binding, loses inhibition" as an unqualified loss.
# Deliberately NOT including "stabilit\w*": it matches the noun phrase
# "subunit-stability measurement", which names a measurement rather than
# asserting anything about the variant.
INTACT = (r"(intact\w*|WT-?like|wild-?type-?like|unaffected|unimpaired|"
          r"no (?:detectable )?(?:effect|change|loss)|normal\w*|preserv\w*|"
          r"retain\w*|unchanged|comparable to (?:WT|wild)|similar to (?:WT|wild)|"
          r"nativ\w*|stable\b)")
LOSS = (r"(abolish\w*|lost|loss|reduc\w*|impair\w*|disrupt\w*|weaken\w*|decreas\w*|"
        r"destabili\w*|eliminat\w*|compromis\w*|defect\w*|diminish\w*)")
# A negated loss ("not interface loss", "no binding defect") is not a loss
# assertion. Stripped before matching so the curator can say what a variant
# does NOT do without tripping the gate.
NEGATED_LOSS = r"\b(?:not|no|without)\s+(?:\w+\s+){0,2}?" + LOSS

# A clause is the unit of assertion. Split on sentence and clause boundaries so
# "WT-like RAF1 binding; impaired GTP hydrolysis" is read as two claims.
CLAUSE_SPLIT = r"[;.]|\s+(?:but|whereas|while|although|though)\s+"


def clauses(basis: str) -> list:
    return [c.strip() for c in re.split(CLAUSE_SPLIT, str(basis)) if c.strip()]


PROXIMITY_WORDS = 4


def _near(text: str, pat: str, vocab: str) -> bool:
    """Does a `pat` match sit within PROXIMITY_WORDS of an axis-vocabulary word?

    Clause-level co-occurrence is too loose: a clause can name the axis and,
    separately, assert something about a different property. Requiring the
    assertion to be grammatically adjacent to the axis noun is what
    distinguishes "WT-like RAF1 binding" (an assertion about binding) from
    "not an independent subunit-stability measurement" (a description of a
    method that happens to contain both).
    """
    toks = re.findall(r"\S+", text)
    hits = [i for i, t in enumerate(toks) if re.search(pat, t, re.I)]
    axes = [i for i, t in enumerate(toks) if re.search(vocab, t, re.I)]
    return any(abs(h - a) <= PROXIMITY_WORDS for h in hits for a in axes)


def assertions_for_axis(basis: str, axis: str):
    """Return (says_intact, says_lost) for assertions attached to this axis."""
    vocab = AXIS_VOCAB.get(axis)
    if vocab is None:
        return False, False
    intact = lost = False
    for c in clauses(basis):
        if not re.search(vocab, c, re.I):
            continue
        if _near(c, INTACT, vocab):
            intact = True
        if _near(re.sub(NEGATED_LOSS, " ", c, flags=re.I), LOSS, vocab):
            lost = True
    return intact, lost


def audit(led: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in led.iterrows():
        tok = str(r.get("expected_token", "")).strip().lower()
        basis = str(r.get("evidence_basis", ""))
        intact, lost = assertions_for_axis(basis, str(r.get("axis", "")))
        conflict = ""
        if tok == "destab" and intact and not lost:
            conflict = "token destab; basis asserts the axis is intact"
        elif tok == "destab" and intact and lost:
            # Both present in axis-relevant clauses: a real contradiction only
            # if the intactness clause is the one carrying the axis vocabulary.
            conflict = "token destab; basis asserts BOTH intact and lost on this axis"
        elif tok == "neutral" and lost and not intact:
            conflict = "token neutral; basis asserts loss on this axis"
        rows.append({"system": r.get("system"), "variant": r.get("variant"),
                     "axis": r.get("axis"), "expected_token": tok,
                     "axis_says_intact": intact, "axis_says_lost": lost,
                     "conflict": conflict, "evidence_basis": basis})
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", default=str(LEDGER))
    ap.add_argument("--verbose", action="store_true",
                    help="print every row, not only conflicts")
    args = ap.parse_args()

    led = pd.read_csv(args.ledger)
    missing = [c for c in ("system", "variant", "axis", "expected_token",
                           "evidence_basis") if c not in led.columns]
    if missing:
        print("FAIL: ledger is missing required columns: %s" % missing)
        return 1

    blank = led["evidence_basis"].isna() | (led["evidence_basis"].astype(str).str.strip() == "")
    if blank.any():
        print("FAIL: %d committed expectation(s) carry no evidence basis" % int(blank.sum()))
        for _, r in led[blank].iterrows():
            print("  - %s %s %s" % (r["system"], r["variant"], r["axis"]))
        return 1

    res = audit(led)
    bad = res[res.conflict != ""]

    if args.verbose:
        for _, r in res.iterrows():
            print("  %-13s %-8s %-13s %-8s intact=%-5s lost=%-5s %s"
                  % (r.system, r.variant, r.axis, r.expected_token,
                     r.axis_says_intact, r.axis_says_lost, r.conflict))

    if len(bad):
        print("FAIL: %d of %d committed expectations contradict their own evidence basis"
              % (len(bad), len(res)))
        for _, r in bad.iterrows():
            print("  - %s %s (%s): %s" % (r.system, r.variant, r.axis, r.conflict))
            print("      token: %s" % r.expected_token)
            print("      basis: %s" % r.evidence_basis[:150])
        print("\nEither the token or the basis is wrong. Do not update the token to "
              "match the prediction; establish which the source supports.")
        return 1

    print("PASS: all %d committed expectations agree with their own evidence basis"
          % len(res))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
