#!/usr/bin/env python
"""Mutation-test the manuscript audits: does each gate actually BIND its value?

An audit that passes on correct prose proves nothing on its own -- it must FAIL
on corrupted prose. This script corrupts one gated value at a time in a scratch
copy of the manuscript and asserts the responsible audit rejects it.

WHY IT EXISTS
=============
Three separate defects of the same species were found by hand in this project:
a check whose expectation was derived from the expression it tested, and twice a
check that could not see prose overstating a value's precision because
``"0.903" in "0.9039"``. When the containment defect was measured
systematically, **45 of 45 digit-ending gates across both audits were
non-binding.** Fixing them by hand is not enough; without this script the
property silently regresses on the next edit.

Two corruption modes, both of which happened for real:

* **overstate precision** -- append a digit to a value where prose states it
  (``0.769`` -> ``0.7699``); the number reported is one no analysis produced.
* **partial update** -- corrupt only ONE mention of a value stated twice. A
  presence-only gate is satisfied by the surviving good mention. This is the
  shape of the drift that left a stale 4,000-draw bound in §3.12 after the
  estimator moved to 200,000 draws.

Exit 0 = every gate binds. Exit 1 = at least one gate is decoration.
"""
from __future__ import annotations

import json
import pathlib
import shutil
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parent.parent
MS_REL = "docs/COMAVI_manuscript_v22.md"
SEC = ("### 3.12 ", "### 3.13 ")
ABSTRACT = ("## Abstract", "## 1. Introduction")


def needles(script: str, anchor: str, var: str, work: pathlib.Path):
    """Recover the literals an audit demands, by dumping them from the audit itself.

    Re-listing them here would create a second source of truth that could drift
    from the audit -- the exact failure this whole file exists to prevent.
    """
    src = (work / "scripts" / f"{script}.py").read_text()
    if src.count(anchor) != 1:
        raise SystemExit(f"{script}: anchor found {src.count(anchor)}x, expected 1 "
                         "(the audit was restructured; update this script)")
    out = work / f"_needles_{script}.json"
    inst = work / "scripts" / f"_inst_{script}.py"
    inst.write_text(src.replace(
        anchor,
        f"    import json as _j, pathlib as _p; "
        f"_p.Path({str(out)!r}).write_text(_j.dumps({var}))\n" + anchor, 1))
    subprocess.run([sys.executable, str(inst)], cwd=work, capture_output=True, text=True)
    inst.unlink()
    if not out.exists():
        raise SystemExit(f"{script}: could not recover gated literals")
    return json.loads(out.read_text())


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        work = pathlib.Path(td) / "repo"
        shutil.copytree(REPO, work, ignore=shutil.ignore_patterns(
            ".git", "_work", "verification_output", "__pycache__", "*.pyc"))
        ms = work / MS_REL
        if not ms.exists():
            print(f"SKIP: {MS_REL} is withheld from this checkout; "
                  "audit-binding cannot be tested here.")
            return 0
        orig = ms.read_text()

        def audit(script: str) -> bool:
            """True if the audit REJECTS the current manuscript."""
            return subprocess.run([sys.executable, f"scripts/{script}.py"], cwd=work,
                                  capture_output=True, text=True).returncode != 0

        specs = [
            ("audit_evidence_claims", "    fails = []", "checks", None),
            ("audit_tier_energy_claims", "    for needle, label in required:", "required", SEC),
        ]

        failures, tested = [], 0
        for script, anchor, var, scope in specs:
            gates = needles(script, anchor, var, work)
            if scope is None:
                lo, hi = 0, len(orig)
            else:
                lo = orig.index(scope[0])
                hi = orig.index(scope[1], lo)
            span = orig[lo:hi]

            for needle, label in gates:
                n = str(needle)
                if not n[-1:].isdigit() or n not in span:
                    continue      # only digit-ending values can be over-precised
                tested += 1
                occ = span.count(n)

                # Mode 1: overstate precision at EVERY mention in scope.
                ms.write_text(orig[:lo] + span.replace(n, n + "9") + orig[hi:])
                if not audit(script):
                    failures.append(f"{script}: gate {label!r} ({n!r}) does not bind -- "
                                    "prose can overstate its precision and the audit passes")

                # Mode 2: corrupt ONE mention only, leaving a correct one behind.
                if occ > 1:
                    k = span.find(n)
                    ms.write_text(orig[:lo] + span[:k + len(n)] + "9" + span[k + len(n):] + orig[hi:])
                    if not audit(script):
                        failures.append(
                            f"{script}: gate {label!r} ({n!r}, stated {occ}x) tolerates a "
                            "partial update -- one mention can go stale undetected")
                ms.write_text(orig)

        # The Abstract restates headline values; a claim reversal there must fail
        # even though every number stays intact.
        a0, a1 = orig.index(ABSTRACT[0]), orig.index(ABSTRACT[1])
        rev = ("none of the", "some of the")
        if rev[0] in orig[a0:a1]:
            tested += 1
            ms.write_text(orig[:a0] + orig[a0:a1].replace(*rev, 1) + orig[a1:])
            if not audit("audit_tier_energy_claims"):
                failures.append("audit_tier_energy_claims: the Abstract's zero-cell claim can "
                                "be REVERSED with every number intact and the audit passes")
            ms.write_text(orig)

        assert ms.read_text() == orig, "scratch manuscript not restored"

        for f in failures:
            print(f"  NON-BINDING  {f}")
        print(f"\naudit-binding: {tested} gates mutation-tested, {len(failures)} non-binding")
        if failures:
            print("FAIL: at least one audit gate is decoration, not a check", file=sys.stderr)
            return 1
        print("PASS: every gated value is bound -- corrupting it fails the audit")
        return 0


if __name__ == "__main__":
    sys.exit(main())
