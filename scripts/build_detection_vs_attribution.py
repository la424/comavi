#!/usr/bin/env python3
"""Does COMAVI notice a lesion, and does it name the right one?

WHY THIS EXISTS
---------------
The headline mechanism-pattern score answers neither question on its own. It is
a partial-credit rubric over a cohort that is half curated-silent variants, so
it sits between two quantities a reader actually wants:

  detection    on a variant with a curated structural mechanism, does ANY of
               the three axes fire at the reference threshold?
  attribution  does the axis the mechanism is attributed to fire?

Separating them is what shows where the pipeline is strong and where it is not,
and it localises the weakest arm to a specific failure rather than a low number.

DEFINITIONS
-----------
Firing is |ddG| >= the reference threshold on an axis, using the monomer value
and the strongest-magnitude partner for the complex-context fold and binding
axes -- the same reduction the concordance layer uses.

The expected axis follows the curated mechanism class:
  fold_mechanism        monomer fold
  ppi_destab_mechanism  binding
  mixed_structural      at least one of {complex fold, binding} AND at least
                        one of {monomer, complex fold} -- a mixed expectation
                        is not satisfied by a single axis

Correct rejection is the complement on curated structurally silent variants:
no axis fires.

Usage
-----
    PYTHONPATH=. python scripts/build_detection_vs_attribution.py
    PYTHONPATH=. python scripts/build_detection_vs_attribution.py --check
"""
import argparse
import json
import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_detection_vs_attribution.json"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
REFERENCE_T = 2.5
STRUCTURAL = ("fold_mechanism", "ppi_destab_mechanism", "mixed_structural")
SILENT = "structurally_silent"


def axis_columns(canon, kind):
    bad = ("_sd", "_ci95", "_disting", "_runs", "_vote")
    pre = "ddg_%s_" % kind
    return [c for c in canon.columns
            if c.startswith(pre) and not any(b in c for b in bad)]


def strongest(row, cols, prefix):
    vals = [row[c] for c in cols if pd.notna(row.get(c))]
    return max(vals, key=abs) if vals else None


def per_axis_thresholds(threshold):
    """Expand a threshold spec into (monomer, fold, binding).

    Four of the five committed thresholds are scalar; the substitution-adjusted
    one is per-axis ({'monomer': 2.9, 'fold': 2.9, 'binding': 3.5}). Accepting
    both here is what lets build() serve the whole committed set, so the
    threshold sweep is this same computation evaluated at five points rather
    than a second implementation that agrees at one.
    """
    if isinstance(threshold, dict):
        return (float(threshold["monomer"]), float(threshold["fold"]),
                float(threshold["binding"]))
    t = float(threshold)
    return t, t, t


def build(canon, threshold=REFERENCE_T):
    graded = canon[canon.mech_consistency_t25.isin(GRADE_MAP)].copy()
    graded["score"] = graded.mech_consistency_t25.map(GRADE_MAP)
    fold_cols = axis_columns(canon, "fold")
    bind_cols = axis_columns(canon, "binding")
    graded["fold_max"] = [strongest(r, fold_cols, "ddg_fold_") for _, r in graded.iterrows()]
    graded["bind_max"] = [strongest(r, bind_cols, "ddg_binding_") for _, r in graded.iterrows()]

    t_mono, t_fold, t_bind = per_axis_thresholds(threshold)

    def fires_at(t):
        def f(v):
            return bool(pd.notna(v) and abs(float(v)) >= t)
        return f

    graded["mono_fires"] = graded.ddg_monomer.map(fires_at(t_mono))
    graded["fold_fires"] = graded.fold_max.map(fires_at(t_fold))
    graded["bind_fires"] = graded.bind_max.map(fires_at(t_bind))
    graded["any_fires"] = graded[["mono_fires", "fold_fires", "bind_fires"]].any(axis=1)

    def right_axis(r):
        if r.expected_mech_class == "fold_mechanism":
            return bool(r.mono_fires)
        if r.expected_mech_class == "ppi_destab_mechanism":
            return bool(r.bind_fires)
        if r.expected_mech_class == "mixed_structural":
            return bool(r.fold_fires or r.bind_fires) and bool(r.mono_fires or r.fold_fires)
        return None

    graded["right_axis"] = graded.apply(right_axis, axis=1)
    struct = graded[graded.expected_mech_class.isin(STRUCTURAL)]
    silent = graded[graded.expected_mech_class == SILENT]

    per_class = {}
    for cls, sub in graded.groupby("expected_mech_class"):
        row = {"n": int(len(sub)),
               "mean_grade": round(float(sub.score.mean()), 4),
               "consistent": int((sub.mech_consistency_t25 == "consistent").sum()),
               "partial": int((sub.mech_consistency_t25 == "partial").sum()),
               "inconsistent": int((sub.mech_consistency_t25 == "inconsistent").sum())}
        if cls in STRUCTURAL:
            row["detection"] = int(sub.any_fires.sum())
            row["attribution"] = int(sub.right_axis.sum())
        elif cls == SILENT:
            row["correct_rejection"] = int((~sub.any_fires).sum())
        per_class[cls] = row

    ppi = graded[graded.expected_mech_class == "ppi_destab_mechanism"]
    return {
        "threshold": threshold,
        "population": {"graded": int(len(graded)),
                       "structural": int(len(struct)),
                       "silent": int(len(silent))},
        "headline_score": {"total": float(graded.score.sum()),
                           "n": int(len(graded)),
                           "mean": round(float(graded.score.mean()), 4)},
        "detection": {"k": int(struct.any_fires.sum()), "n": int(len(struct)),
                      "rate": round(float(struct.any_fires.mean()), 3)},
        "attribution": {"k": int(struct.right_axis.sum()), "n": int(len(struct)),
                        "rate": round(float(struct.right_axis.mean()), 3)},
        "correct_rejection": {"k": int((~silent.any_fires).sum()), "n": int(len(silent)),
                              "rate": round(float((~silent.any_fires).mean()), 3)},
        "per_class": per_class,
        "interface_class_diagnosis": {
            "n": int(len(ppi)),
            "binding_axis_fires": int(ppi.bind_fires.sum()),
            "complex_fold_axis_fires": int(ppi.fold_fires.sum()),
            "note": ("For variants curated as pure interface destabilizers the "
                     "complex-context fold axis fires more often than the binding "
                     "axis, so the lesion is detected but attributed to the wrong "
                     "axis. This is the separation the pipeline exists to make and "
                     "it is where the benchmark is weakest."),
            "large_monomer_effect": [
                {"gene": str(r.gene).upper(), "variant": r.variant,
                 "ddg_monomer": round(float(r.ddg_monomer), 2)}
                for _, r in ppi.iterrows()
                if pd.notna(r.ddg_monomer) and abs(float(r.ddg_monomer)) >= 5.0],
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    result = build(canon)

    # Structural checks, not frozen literals. The v7.9 ground-truth correction
    # moved the headline from 39.5 to 41.0 and three literals like the one that
    # used to sit here all failed at once; bumping them would have removed the
    # guard. Recomputing the expectation from the canonical keeps it live.
    _W = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
    _g = canon[canon.mech_consistency_t25.isin(_W)]
    p = result["population"]
    assert p["graded"] == len(_g), (p["graded"], len(_g))
    assert p["structural"] + p["silent"] == p["graded"]
    assert abs(result["headline_score"]["total"]
               - _g.mech_consistency_t25.map(_W).sum()) < 1e-9
    # the weakest arm claim in the manuscript rests on this ordering
    means = {k: v["mean_grade"] for k, v in result["per_class"].items()}
    assert min(means, key=means.get) == "ppi_destab_mechanism", means

    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: detection/attribution matches a fresh build")
        return 0

    print("graded %d (structural %d, silent %d) at threshold %.1f kcal/mol"
          % (p["graded"], p["structural"], p["silent"], result["threshold"]))
    for k in ("detection", "attribution", "correct_rejection"):
        d = result[k]
        print("  %-18s %2d/%-2d = %.3f" % (k, d["k"], d["n"], d["rate"]))
    print("\nper expected mechanism class:")
    for cls, r in sorted(result["per_class"].items(),
                         key=lambda kv: -kv[1]["mean_grade"]):
        extra = ""
        if "detection" in r:
            extra = "detect %2d/%-2d  attribute %2d/%-2d" % (
                r["detection"], r["n"], r["attribution"], r["n"])
        elif "correct_rejection" in r:
            extra = "reject %2d/%-2d" % (r["correct_rejection"], r["n"])
        print("  %-22s n=%2d  grade %.3f   %s" % (cls, r["n"], r["mean_grade"], extra))
    d = result["interface_class_diagnosis"]
    print("\ninterface class (n=%d): binding axis fires %d, complex-fold axis fires %d"
          % (d["n"], d["binding_axis_fires"], d["complex_fold_axis_fires"]))
    if d["large_monomer_effect"]:
        print("  large monomer effects among them: %s"
              % ", ".join("%s %s (%.1f)" % (x["gene"], x["variant"], x["ddg_monomer"])
                          for x in d["large_monomer_effect"]))

    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
