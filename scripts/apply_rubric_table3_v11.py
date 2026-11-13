#!/usr/bin/env python3
"""Re-grade the canonical with the grader that implements Table 3 as published.

WHY THIS EXISTS
---------------
v11 corrected two places where grade_mechanism_consistency (apply_concordance_v5)
departed from Table 3, the rubric the manuscript publishes:

  1. Partial credit without recovering anything. Table 3 grades an expected
     effect missed entirely as inconsistent; the code gave 0.5 to every miss
     without a false call on a tested axis, including calls asserting no
     energetic effect and calls made only on axes the literature never
     committed.
  2. No direction check for mixed-direction labels. Table 3 grades a predicted
     change opposite in direction to the expected effect as inconsistent; the
     code skipped the check for every "(conflicting)" label.

The canonical's grade columns were written by the old grader, so they are
re-derived here from the stored mechanism calls: for both pipelines (mech_*,
nbhd_mech_*) at all five thresholds, with the false-positive and
missed-positive lists, the t=2.5 summary and the threshold-stability flag.
Only cells whose value changes are written; every other byte of the table is
left as it was, and the run asserts that.

--check re-derives every grade with the current grader and fails if the
canonical disagrees, so the stored grades can no longer drift from the code.

Usage:
    python scripts/apply_rubric_table3_v11.py           # apply, write the audit
    python scripts/apply_rubric_table3_v11.py --check   # verify only
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
AUDIT = REPO / "docs" / "audit" / "COMAVI_v11_rubric_regrade.csv"
TAGS = ("t10", "t15", "t20", "t25", "tSAP")
PIPELINES = (("", "mech_"), ("nbhd_", "nbhd_mech_"))
GRADES = ("consistent", "partial", "inconsistent")


def blank(v):
    return v is None or (isinstance(v, float) and pd.isna(v)) or v == ""


def as_cell(v):
    """The canonical stores 'N/A' and empty lists as empty cells."""
    return float("nan") if (v in ("", "N/A") or v is None) else v


def same(a, b):
    return (blank(a) and blank(b)) or a == b


def reason(row, label, old, new, fp, mp):
    if new == "inconsistent" and old == "partial":
        fold_m, fold_c, bind, _i, _d = ac.MECH_LABEL_ASSERTIONS.get(label, (0, 0, 0, 0, ""))
        if not (fold_m or fold_c or bind):
            return "expected effect missed entirely; the call asserts no energetic effect"
        return "expected effect missed entirely; the call names only axes the literature does not commit"
    if new == "inconsistent" and old == "consistent":
        return "predicted direction opposite to the expected effect on an asserted axis"
    return "other"


def regrade(canon: pd.DataFrame):
    out = canon.copy()
    changes = []
    for idx, row in canon.iterrows():
        cls = row.get("expected_mech_class")
        if blank(cls):
            continue
        # Rows excluded from grading by a recorded decision (BRCA1 R1699L/Q,
        # mechanism_not_observable_on_structure) stay ungraded: the exclusion
        # was applied after grading, so the grader alone would grade them.
        if not blank(row.get("mech_graded_excluded_reason")):
            continue
        axes = ac.classify_axis_status(row)
        if canon.at[idx, "mech_consistency_t25"] in GRADES:
            derived = ac.derive_expected_mech_class(row)
            if derived != cls:
                raise SystemExit("STOP: %s %s stored class %r disagrees with derived %r"
                                 % (row.system, row.variant, cls, derived))
        for prefix, call_prefix in PIPELINES:
            for tag in TAGS:
                label = row.get("%s%s" % (call_prefix, tag))
                g, fp, mp = ac.grade_mechanism_consistency(row, label, cls, axes)
                cells = {
                    "%smech_consistency_%s" % (prefix, tag): as_cell(g),
                    "%smech_false_positive_axes_%s" % (prefix, tag): as_cell(",".join(fp) if fp else ""),
                    "%smech_missed_positive_axes_%s" % (prefix, tag): as_cell(",".join(mp) if mp else ""),
                }
                old_grade = canon.at[idx, "%smech_consistency_%s" % (prefix, tag)]
                for col, val in cells.items():
                    if not same(canon.at[idx, col], val):
                        out.at[idx, col] = val
                if not same(old_grade, as_cell(g)):
                    changes.append(dict(
                        pipeline="primary" if not prefix else "neighbourhood",
                        threshold=tag, system=row.system, variant=row.variant,
                        expected_class=cls, call=label, old_grade=old_grade, new_grade=g,
                        missed_axes=",".join(mp), reason=reason(row, label, old_grade, g, fp, mp)))
        for prefix, _ in PIPELINES:
            grades = [out.at[idx, "%smech_consistency_%s" % (prefix, t)] for t in TAGS]
            if any(gr in GRADES for gr in grades):
                summ = out.at[idx, "%smech_consistency_t25" % prefix]
                if not same(out.at[idx, "%smech_consistency_summary" % prefix], summ):
                    out.at[idx, "%smech_consistency_summary" % prefix] = summ
                stable = len(set(grades)) == 1
                col = "%smech_consistency_threshold_stable" % prefix
                if not same(out.at[idx, col], stable):
                    out.at[idx, col] = stable
    return out, pd.DataFrame(changes)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    canon = pd.read_csv(CANON, low_memory=False)
    new, changes = regrade(canon)
    if args.check:
        if len(changes):
            print(changes[["pipeline", "threshold", "variant", "old_grade", "new_grade"]].to_string(index=False))
            print("FAIL: %d stored grades disagree with the Table 3 grader" % len(changes))
            return 1
        print("PASS: every stored grade, both pipelines, all five thresholds, matches the Table 3 grader")
        return 0
    touched = [c for c in canon.columns if not canon[c].equals(new[c])]
    allowed = tuple("%smech_%s" % (p, s) for p, _ in PIPELINES
                    for s in ("consistency_", "false_positive_axes_", "missed_positive_axes_"))
    stray = [c for c in touched if not c.startswith(allowed)]
    if stray:
        raise SystemExit("STOP: columns outside the grade family would change: %s" % stray)
    new.to_csv(CANON, index=False)
    changes.to_csv(AUDIT, index=False)
    print("re-graded: %d grade cells changed in %d columns" % (len(changes), len(touched)))
    if len(changes):
        print(changes.groupby(["pipeline", "threshold"]).size().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
