#!/usr/bin/env python3
"""Summary records for the two priority-score sensitivity analyses in S2 Text.

WHY THIS EXISTS
---------------
S2 Text makes two claims that nothing in this repository recomputed:

  "Eight of the nine keep a graded context term, and all eight rank the cohort
   almost identically to the reported score (Spearman at least 0.938); their
   ROC AUC spans 0.794 to 0.844"

  "In all 13 folds the combined score ranked above both of its own components"

Both were computed in an interactive session and survive only as saved
tables. The underlying quantities are committed -- analyze_isds_v1.py already
writes all nine alternative scores as columns of ISDS_v1_per_variant.csv, and
the leave-one-system-out folds are a regrouping of that same file -- so what
was missing is only the summary each sentence quotes. A published number that
nothing recomputes is how the repository has acquired stale values before, and
these are the last two in the supplement without a generator.

Nothing here redefines a score. Every alternative formula is read from the
column the pipeline already produced, which is why this file is a summary
rather than a reimplementation: a second implementation could drift from the
pipeline and agree with the paper, which is the failure mode to avoid.

OUTPUT
------
reference_outputs/COMAVI_priority_sensitivity_summary.json
"""
import argparse
import json
import pathlib
import sys

import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

REPO = pathlib.Path(__file__).resolve().parents[1]
PERVAR = REPO / "reference_outputs" / "isds_v1" / "ISDS_v1_per_variant.csv"
OUT = REPO / "reference_outputs" / "COMAVI_priority_sensitivity_summary.json"

# Label -> the column analyze_isds_v1.py writes for that alternative.
ALTERNATIVES = {
    "hard cap": "sens_hard_cap",
    "exponential": "sens_exponential",
    "smooth normalization": "sens_smooth_norm",
    "common per-axis anchor": "sens_common_2_5",
    "energy weight 0.4": "sens_energy_weight_0_4",
    "energy weight 0.6": "sens_energy_weight_0_6",
    "tier without interface bonus": "sens_no_interface_tier",
    "interface-only context": "sens_interface_only_context",
    "cohort-standardized": "historical_cohort_z_combination",
}
# The one alternative that replaces the graded tier with a binary flag. S2
# separates it from the other eight precisely because it is not a graded
# context term, so the "all eight" claim excludes it by definition, not by
# selection after seeing the result.
BINARY_CONTEXT = "interface-only context"
REPORTED = "isds_v1"
COMPONENTS = ("isds_energy_component", "isds_context_component")


def build():
    pv = pd.read_csv(PERVAR)
    missing = [c for c in list(ALTERNATIVES.values()) + [REPORTED] + list(COMPONENTS)
               if c not in pv.columns]
    assert not missing, "per-variant record is missing columns: %s" % missing
    y = pv.structural_ground_truth.astype(bool)

    alts = []
    for label, col in ALTERNATIVES.items():
        alts.append({
            "alternative": label,
            "column": col,
            "graded_context_term": label != BINARY_CONTEXT,
            "roc_auc": round(float(roc_auc_score(y, pv[col])), 4),
            "spearman_vs_reported": round(float(stats.spearmanr(pv[col],
                                                                pv[REPORTED]).statistic), 4),
        })
    alts.sort(key=lambda r: r["roc_auc"])
    graded = [a for a in alts if a["graded_context_term"]]

    folds = []
    for dropped in sorted(pv.system.unique()):
        d = pv[pv.system != dropped]
        yy = d.structural_ground_truth.astype(bool)
        if yy.sum() < 2 or (~yy).sum() < 2:
            continue
        row = {"dropped_system": dropped,
               "combined": round(float(roc_auc_score(yy, d[REPORTED])), 4)}
        for c in COMPONENTS:
            row[c] = round(float(roc_auc_score(yy, d[c])), 4)
        row["combined_above_both"] = bool(row["combined"] > row[COMPONENTS[0]]
                                          and row["combined"] > row[COMPONENTS[1]])
        folds.append(row)

    out = {
        "population": {"n": int(len(pv)), "positive": int(y.sum()),
                       "negative": int((~y).sum()),
                       "systems": int(pv.system.nunique())},
        "reported_roc_auc": round(float(roc_auc_score(y, pv[REPORTED])), 4),
        "alternative_formulas": {
            "n": len(alts),
            "rows": alts,
            "graded_context_n": len(graded),
            "graded_context_auc_span": [min(a["roc_auc"] for a in graded),
                                        max(a["roc_auc"] for a in graded)],
            "graded_context_min_spearman": min(a["spearman_vs_reported"] for a in graded),
            "binary_context": next(a for a in alts if not a["graded_context_term"]),
            "note": ("All nine are read from the columns the pipeline writes; "
                     "none is reimplemented here."),
        },
        "leave_one_system_out": {
            "n_folds": len(folds),
            "folds": folds,
            "combined_above_both_in": sum(f["combined_above_both"] for f in folds),
        },
        "caveat": ("Both analyses are in-sample: the same benchmark informed "
                   "the development of the reported score."),
    }

    # The two S2 sentences depend on these shapes, not merely on these values.
    los = out["leave_one_system_out"]
    assert los["combined_above_both_in"] == los["n_folds"], (
        "the combined score no longer beats both components in every fold; "
        "the S2 sentence must be rewritten")
    assert out["alternative_formulas"]["graded_context_n"] == 8, (
        "the count of graded-context alternatives changed; S2 says eight")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    fresh = json.loads(json.dumps(build()))
    if args.check:
        if not OUT.exists():
            print("FAIL: %s is missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != fresh:
            print("FAIL: priority sensitivity summary is stale")
            return 1
        af = fresh["alternative_formulas"]
        print("PASS: priority sensitivity summary reproduces "
              "(graded-context AUC span %s, combined above both in %d of %d folds)"
              % (af["graded_context_auc_span"],
                 fresh["leave_one_system_out"]["combined_above_both_in"],
                 fresh["leave_one_system_out"]["n_folds"]))
        return 0
    OUT.write_text(json.dumps(fresh, indent=1) + "\n")
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
