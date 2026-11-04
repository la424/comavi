#!/usr/bin/env python3
"""Assert the manuscript's predicted-versus-measured claims match their records.

WHY THIS EXISTS
---------------
The calibration family -- the rank correlation, the fitted slope, the image of
2.5 kcal/mol, the Glu73 sensitivity, the SKEMPI subset, the call relationships
and the two within-benchmark anchors -- was the only group of published values
with no audit over the prose. Two of them went stale and nothing caught either:

  * "measured and predicted calls agreed for 43 of the 63 comparisons" while the
    record read 42. One comparison had moved from agreement to substantive
    disagreement when the hemoglobin measured values were harmonized to a
    per-tetramer convention, and the figure caption was updated while the
    Results sentence was not.
  * "hemoglobin assembly energies (rho = 0.90, n = 5)" while the calibration
    points give 0.80. The same harmonization moved HBB N102T from 1.5 to 3.0
    kcal/mol, which breaks a rank pair.

Both are the same shape as every other defect this repository gates against: a
published number whose generator moved underneath it. Every needle below is
BUILT from a committed record and asserted to appear in the extracted
manuscript text, so a record that moves without the prose following it fails
here rather than in review.

Correlations are recomputed from COMAVI_delta_calibration_points.csv rather than
read from a summary, because the two within-benchmark anchors (BRCT and the
hemoglobin tetramer) have no summary field of their own -- which is precisely
why the hemoglobin value could drift unnoticed.

Usage:
    PYTHONPATH=. python scripts/audit_calibration_claims.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import pandas as pd
from scipy import stats

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"
MS = REPO / "docs" / "COMAVI_manuscript_current.txt"


def needles():
    """[(needle, label)] -- every value computed, none typed."""
    stats_rec = json.loads((RO / "COMAVI_delta_calibration_stats.json").read_text())
    eff = json.loads((RO / "COMAVI_measured_effect_tables.json").read_text())
    points = pd.read_csv(RO / "COMAVI_delta_calibration_points.csv")
    g73 = stats_rec["sensitivity_excluding_glu73"]
    classes = {c["call_relationship"]: c["n"] for c in eff["call_relationships"]["classes"]}
    recovered = {r["threshold"]: r["cell"] for r in eff["measured_effects_recovered"]["by_threshold"]}
    pop = eff["measured_effects_recovered"]["population"]

    def rho_of(system):
        d = points[points["system"] == system]
        return stats.spearmanr(d["measured_kcal"], d["foldx_ddg"]).statistic, len(d)

    brct_rho, brct_n = rho_of("BRCA1 BRCT")
    hb_rho, hb_n = rho_of("Hb tetramer")

    out = [
        # ---- pooled calibration
        ("Spearman \u03c1 = %.2f" % stats_rec["rho"], "pooled rank correlation"),
        ("the fitted slope was %.2f" % stats_rec["slope"], "fitted slope"),
        ("predicted %.2f kcal/mol" % stats_rec["image_of_2p5"], "image of 2.5 kcal/mol"),
        ("\u03c1 = \u2212%.2f" % abs(stats_rec["rho_delta_vs_measured"]), "error against measured"),
        ("%d comparisons from four systems" % stats_rec["n_fit"], "pooled n"),
        # ---- Glu73 sensitivity
        ("rises to %.2f (n = %d)" % (g73["rho"], g73["n_fit"]), "Glu73-excluded correlation"),
        ("the slope barely changes (%.2f)" % g73["slope"], "Glu73-excluded slope"),
        # ---- call relationships at the reference threshold
        ("agreed for %d of the %d comparisons"
         % (classes["Agreement"], sum(classes.values())), "call agreement"),
        ("borderline, n = %d" % classes["Borderline disagreement"], "borderline class"),
        ("substantive, n = %d" % classes["Substantive disagreement"], "substantive class"),
        # ---- measured-effect recovery across thresholds
        ("It did for %s at 1.0 kcal/mol, %s at 2.5 kcal/mol and %s at the prediction-interval limit"
         % (recovered["1.0"].split("/")[0], recovered["2.5"].split("/")[0],
            recovered["tSAP"].split("/")[0]), "recovery ladder"),
        ("the %d comparisons with a measured destabilization of at least 1.0 kcal/mol, "
         "%d from benchmark variants and %d from external measurements"
         % (pop["n"], pop["in_benchmark"], pop["external"]), "recovery population"),
        # ---- the two within-benchmark anchors, recomputed from the points
        # The prose says "BRCT variants"; the record's system label is
        # "BRCA1 BRCT". The needle follows the paper's wording -- only the
        # NUMBERS come from the record.
        ("across the %d graded BRCT variants (\u03c1 = %.2f" % (brct_n, brct_rho),
         "BRCT anchor"),
        ("(\u03c1 = %.2f, n = %d; Fig 5D)" % (hb_rho, hb_n), "hemoglobin anchor"),
    ]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--allow-missing-manuscript", action="store_true",
        help="Skip instead of failing when the manuscript text is absent. For "
             "manuscript-withheld clones ONLY; a silent skip here is how the "
             "prose drifted from the data in the first place.")
    args = ap.parse_args()

    if not MS.exists():
        if args.allow_missing_manuscript:
            print("[SKIP] %s absent; skip explicitly requested." % MS.relative_to(REPO))
            return 0
        print("[FAIL] %s is not present. Regenerate it with "
              "scripts/extract_manuscript_text.py." % MS.relative_to(REPO), file=sys.stderr)
        return 1

    text = MS.read_text()
    missing = [(n, lab) for n, lab in needles() if n not in text]
    for n, lab in missing:
        print("  MISSING [%s]: %r" % (lab, n))
    checked = len(needles())
    if missing:
        print("FAIL: %d of %d calibration claims disagree with their records"
              % (len(missing), checked))
        return 1
    print("PASS: all %d calibration claims match their records" % checked)
    return 0


if __name__ == "__main__":
    sys.exit(main())
