#!/usr/bin/env python3
"""Assert the current-column values in RELEASE_NOTES_*.md match the records.

WHY THIS EXISTS
---------------
Release notes must restate headline numbers -- that is what they are for -- so
verify_no_duplicated_headlines.py cannot simply forbid them. But an unchecked
document that restates a headline is exactly how this project has produced drift
before: a figures README quoting a cohort that no longer existed, a results
summary whose template hardcoded a fraction beside the computed rate.

So the notes are allowed there and gated here instead. Every needle below is
BUILT from a committed record and asserted to appear in the notes, in the bold
form the current column uses. If a record moves and the notes are not updated,
this fails. The previous-release column is deliberately not checked against the
records -- those values are superseded by construction and are verified against
what the previous tag actually published, which git holds and these notes
document.

Usage:
    python scripts/verify_release_notes_values.py
    python scripts/verify_release_notes_values.py --notes RELEASE_NOTES_v3.0-methods.md
"""
import argparse
import glob
import json
import pathlib
import re
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"


def expected_current():
    """[(label, needle)] -- every needle computed, none typed."""
    sweep = {r["threshold_tag"]: r for r in
             json.loads((RO / "COMAVI_threshold_sweep.json").read_text())["rows"]}
    t25 = sweep["t25"]
    strat = json.loads((RO / "COMAVI_evidence_stratified_agreement.json").read_text())
    stress = json.loads((RO / "stress_tests" /
                         "COMAVI_stress_verified_statistics.json").read_text())
    led = pd.read_csv(RO / "COMAVI_evidence_ledger.csv")
    pm = pd.read_csv(RO / "isds_v1" / "ISDS_v1_primary_metrics.csv").set_index("score")
    pv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
    can = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
    ranked = (pv.sort_values("isds_v1", ascending=False)
                ["structural_ground_truth"].astype(bool).values)
    npos = int(ranked.sum())
    ci = stress["system_cluster_ci95"]["mc"]["reportable_3"]
    c, acc, st = strat["crude"], strat["accounting"], strat["stratified"]
    e1 = strat["by_evidence_type"]["E1_quantitative_energetic"]

    # (label, needle, in_comparison_table). Table rows must appear in the bold
    # current column; a value present only unbolded there is the superseded one.
    return [
        ("mechanism-pattern agreement",
         "%.4f (%.1f/%d)" % (t25["whole_variant"], t25["whole_variant_k"],
                             t25["whole_variant_n"])),
        ("graded for whole-variant localization",
         "%d of %d" % (t25["whole_variant_n"], len(can))),
        ("committed axis expectations", "%d rows" % len(led)),
        ("four-output continuity aggregate",
         "%d/%d = %.3f" % (t25["structural_agreement_k"],
                           t25["structural_agreement_n"],
                           t25["structural_agreement"])),
        ("prioritization population",
         "%d (%d / %d)" % (len(pv), npos, len(pv) - npos)),
        ("priority-score ROC AUC", "%.3f" % pm.loc["isds_v1", "roc_auc"]),
        ("top-20 recovery", "%d/%d" % (int(ranked[:20].sum()), npos)),
        # Stated in prose rather than in the comparison table, so presence is
        # sufficient: there is no superseded twin on the same row for a bold
        # marker to distinguish it from.
        ("system-cluster 95% interval", "%.3f\u2013%.3f" % (ci[0], ci[1]), False),
        # Supplementary evidence stratification
        ("committed axes / total gradeable",
         "%d / %d" % (acc["n_committed_axes"], acc["n_total_gradeable"])),
        ("E1\u2013E2 quantitative evidence",
         "%d/%d = %.3f" % (c["quant_ok"], c["quant_n"], c["quant_rate"])),
        ("E3\u2013E5 qualitative or inferred",
         "%d/%d = %.3f" % (c["soft_ok"], c["soft_n"], c["soft_rate"])),
        ("crude contrast, Fisher p", "%.4f" % c["p"]),
        ("Mantel-Haenszel OR",
         "%.3f (p = %.3f)" % (st["OR_mantel_haenszel"], st["p"])),
        ("E1 quantitative energetic axes",
         "%d/%d = %.3f" % (e1["ok"], e1["n"], e1["rate"])),
    ]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--notes", default=None,
                    help="a single notes file; default is every RELEASE_NOTES_*.md")
    args = ap.parse_args()
    files = ([REPO / args.notes] if args.notes
             else sorted(pathlib.Path(p) for p in glob.glob(str(REPO / "RELEASE_NOTES_*.md"))))
    if not files:
        print("FAIL: no RELEASE_NOTES_*.md found; this gate is wired into CI and "
              "a zero-file glob must not report success")
        return 1

    needles = expected_current()
    bad = []
    for path in files:
        text = path.read_text()
        for label, needle, *rest in needles:
            in_table = rest[0] if rest else True
            # The current column is bold; a value present only unbolded is in the
            # superseded column and does not satisfy the claim.
            if in_table:
                ok = ("**%s**" % needle) in text or ("**%s" % needle) in text
                why = ("present but not in the bold current column"
                       if needle in text else "absent")
            else:
                ok = needle in text
                why = "absent"
            if not ok:
                bad.append((path.name, label, needle, why))
    for name, label, needle, why in bad:
        print("  MISSING [%s] %s: expected %r, %s" % (name, label, needle, why))
    if bad:
        print("FAIL: %d current-column value(s) disagree with the records (%d checked)"
              % (len(bad), len(needles) * len(files)))
        return 1
    print("PASS: every current-column value in %d release-notes file(s) matches the "
          "records (%d checked)" % (len(files), len(needles) * len(files)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
