#!/usr/bin/env python3
"""What would a trivial strategy score on the mechanism-pattern rubric?

WHY THIS EXISTS
---------------
The manuscript reports a whole-variant mechanism-pattern score of 0.693 without
a reference class, so a reader cannot tell whether that is good. The obvious
worry -- that the number is carried by the curated structurally silent half of
the cohort -- is answerable and the answer is no, but only if the do-nothing
strategy is actually scored rather than argued about.

The cohort is near-balanced (29 curated silent, 28 committed structural), so
0.693 is already a balanced measure: the unweighted mean of the two arms is
0.691, within 0.003 of the weighted mean. What it needs is not a caveat but a
floor.

METHOD
------
Each trivial strategy emits one fixed mechanism label for every variant, and
that label is graded by the SHIPPED rubric -- apply_concordance_v5
.grade_mechanism_consistency, with axis status from classify_axis_status -- on
the same 57 graded variants and the same expected classes COMAVI is scored
against. Nothing is reimplemented, so a rubric change moves the baselines and
the observed score together.

The strongest trivial strategy is "never fire", which scores every curated
silent variant correct and every structural variant wrong. It is the number to
beat and the honest floor to quote.

Usage
-----
    PYTHONPATH=.:scripts python scripts/build_rubric_baselines.py
    PYTHONPATH=.:scripts python scripts/build_rubric_baselines.py --check
"""
import argparse
import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_rubric_baselines.json"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
STRATEGIES = [
    ("never_fire", "No structural effect detected",
     "Calls every variant structurally silent. Scores every curated silent "
     "variant correct and misses every structural mechanism."),
    ("always_fold_and_ppi", "Both fold + PPI destabilization",
     "Calls every variant destabilizing on both fold subtypes and binding."),
    ("always_monomer_fold", "Monomer fold destabilization",
     "Calls every variant a monomer-fold destabilizer."),
    ("always_multimer_fold", "Multimer fold destabilization",
     "Calls every variant a complex-context fold destabilizer."),
    ("always_ppi", "PPI destabilization",
     "Calls every variant an interface destabilizer."),
]


def build(canon):
    graded = canon[canon.mech_consistency_t25.isin(GRADE_MAP)].copy()
    axes = {i: ac.classify_axis_status(r) for i, r in graded.iterrows()}

    def score_fixed_label(label):
        vals, grades = [], []
        for i, r in graded.iterrows():
            g, _, _ = ac.grade_mechanism_consistency(
                r, label, r.get("expected_mech_class"), axes[i])
            grades.append(g)
            if g in GRADE_MAP:
                vals.append(GRADE_MAP[g])
        return vals, grades

    baselines = {}
    for key, label, note in STRATEGIES:
        vals, grades = score_fixed_label(label)
        baselines[key] = {
            "mechanism_label": label,
            "n_graded": len(vals),
            "score": round(sum(vals) / len(vals), 4),
            "grade_counts": {g: grades.count(g) for g in sorted(set(grades))},
            "note": note,
        }

    observed = float(graded.mech_consistency_t25.map(GRADE_MAP).mean())
    silent = graded.expected_mech_class == "structurally_silent"
    arm_silent = float(graded[silent].mech_consistency_t25.map(GRADE_MAP).mean())
    arm_struct = float(graded[~silent].mech_consistency_t25.map(GRADE_MAP).mean())
    best_key = max(baselines, key=lambda k: baselines[k]["score"])

    return {
        "population": {"n_graded": int(len(graded)),
                       "curated_silent": int(silent.sum()),
                       "committed_structural": int((~silent).sum()),
                       "silent_fraction": round(float(silent.mean()), 4)},
        "observed": {"score": round(observed, 4),
                     "arm_curated_silent": round(arm_silent, 4),
                     "arm_committed_structural": round(arm_struct, 4),
                     "unweighted_mean_of_arms": round((arm_silent + arm_struct) / 2, 4),
                     "note": ("The cohort is near-balanced, so the weighted score and the "
                              "unweighted mean of the two arms agree to within 0.003. The "
                              "headline is already a balanced measure, not an average "
                              "dominated by the easier class.")},
        "baselines": baselines,
        "best_trivial_strategy": {
            "key": best_key,
            "score": baselines[best_key]["score"],
            "margin_over_best_trivial": round(observed - baselines[best_key]["score"], 4),
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    result = build(canon)

    p, o = result["population"], result["observed"]
    assert p["n_graded"] == 57, p["n_graded"]
    assert abs(o["score"] - o["unweighted_mean_of_arms"]) < 0.01, o
    # the claim the manuscript will make
    assert result["best_trivial_strategy"]["key"] == "never_fire", result["best_trivial_strategy"]
    assert result["best_trivial_strategy"]["margin_over_best_trivial"] > 0

    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: rubric baselines match a fresh build")
        return 0

    print("population: %d graded (%d curated silent, %d committed structural)"
          % (p["n_graded"], p["curated_silent"], p["committed_structural"]))
    print("\nscores on the shipped rubric, ascending:")
    rows = sorted(((v["score"], k) for k, v in result["baselines"].items()))
    for s, k in rows:
        print("  %-24s %.4f  %s" % (k, s, "#" * int(s * 44)))
    print("  %-24s %.4f  %s" % ("COMAVI (observed)", o["score"], "#" * int(o["score"] * 44)))
    b = result["best_trivial_strategy"]
    print("\nbest trivial strategy is %s at %.4f; COMAVI exceeds it by %.4f"
          % (b["key"], b["score"], b["margin_over_best_trivial"]))
    print("arms: curated silent %.4f, committed structural %.4f (unweighted mean %.4f)"
          % (o["arm_curated_silent"], o["arm_committed_structural"],
             o["unweighted_mean_of_arms"]))

    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
