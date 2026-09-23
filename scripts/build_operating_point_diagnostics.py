#!/usr/bin/env python3
"""Operating-point diagnostics: class mixture, reweighting, paired resampling.

WHY THIS EXISTS
---------------
S4 Text reports three quantities about how the choice of reference threshold
depends on the cohort rather than on the method:

  1. the benchmark's structural fraction (class mixture);
  2. the structural-class weight at which 2.5 kcal/mol stops being the pooled
     optimum, which says how much of the apparent optimum is class mixture;
  3. a paired whole-system bootstrap of graded recovery and correct rejection
     at each threshold, with the probability that each moves in the stated
     direction.

None of them had a generator. The consequence was that the whole block went
stale at the v7.9 ground-truth correction and survived an audit that checks
every decimal against the corpus of generator-emitted values, because each
stale value happens to occur somewhere in a large committed CSV. All four
paired-resample endpoints reproduced exactly from the pre-correction 28/29
class split, and the block stated the structural fraction twice with two
different values -- 0.439 (25/57, correct) and 0.491 (28/57, stale) -- three
sentences apart.

Paired means the same system resample is scored at every threshold, so the
threshold contrast is within-resample and the probability of a direction is
not inflated by between-resample variance in the cohort composition.

Usage:
    PYTHONPATH=. python scripts/build_operating_point_diagnostics.py
    PYTHONPATH=. python scripts/build_operating_point_diagnostics.py --check
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_operating_point_diagnostics.json"

GRADE = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
SILENT = "structurally_silent"
TAGS = ["t10", "t15", "t20", "t25", "tSAP"]
KCAL = {"t10": "1.0", "t15": "1.5", "t20": "2.0", "t25": "2.5",
        "tSAP": "2.9/2.9/3.5"}
N_BOOT = 10000
SEED = 11


def graded(canon):
    cols = ["mech_consistency_%s" % tg for tg in TAGS]
    missing = [c for c in cols if c not in canon.columns]
    if missing:
        raise SystemExit("canonical lacks graded columns: %s" % missing)
    g = canon[canon.mech_consistency_t25.isin(GRADE)].copy()
    for tg, c in zip(TAGS, cols):
        g["grade_%s" % tg] = g[c].map(GRADE)
    g["is_silent"] = g.expected_mech_class.eq(SILENT)
    return g


def arms(sub):
    """Graded recovery among structural rows, correct rejection among silent."""
    out = {}
    s = sub[~sub.is_silent]
    q = sub[sub.is_silent]
    for tg in TAGS:
        out[tg] = (float(s["grade_%s" % tg].mean()) if len(s) else float("nan"),
                   float(q["grade_%s" % tg].mean()) if len(q) else float("nan"))
    return out


def build(canon):
    g = graded(canon)
    n_struct = int((~g.is_silent).sum())
    n_silent = int(g.is_silent.sum())
    frac = n_struct / len(g)

    point = arms(g)

    # [2] class reweighting. Pooled score under weight w on the structural arm.
    # Report the largest w at which t25 is still the argmax, on a fine grid.
    grid = np.round(np.arange(0.0, 1.0001, 0.001), 3)
    still = []
    for w in grid:
        pooled = {tg: w * point[tg][0] + (1 - w) * point[tg][1] for tg in TAGS}
        if max(pooled, key=pooled.get) == "t25":
            still.append(float(w))
    w_hi = max(still) if still else None
    w_lo = min(still) if still else None

    # [3] paired whole-system bootstrap: resample systems, score every
    # threshold on the SAME resample.
    rng = np.random.default_rng(SEED)
    systems = g.system.unique()
    idx = {s: g.index[g.system.eq(s)].to_numpy() for s in systems}
    rec = {tg: np.empty(N_BOOT) for tg in TAGS}
    rej = {tg: np.empty(N_BOOT) for tg in TAGS}
    for b in range(N_BOOT):
        pick = rng.choice(systems, size=len(systems), replace=True)
        rows = np.concatenate([idx[s] for s in pick])
        a = arms(g.loc[rows])
        for tg in TAGS:
            rec[tg][b], rej[tg][b] = a[tg]

    def ci(v):
        v = v[~np.isnan(v)]
        return [round(float(np.percentile(v, 2.5)), 3),
                round(float(np.percentile(v, 97.5)), 3)]

    lo, hi = TAGS[0], TAGS[-1]
    p_rec_fell = float(np.mean(rec[hi] < rec[lo]))
    p_rej_rose = float(np.mean(rej[hi] > rej[lo]))

    return {
        "population": {"n_graded": len(g), "n_structural": n_struct,
                       "n_silent": n_silent,
                       "structural_fraction": round(frac, 4),
                       "n_systems": int(len(systems))},
        "point_estimates": {tg: {"threshold_kcal_mol": KCAL[tg],
                                 "recovery": round(point[tg][0], 4),
                                 "correct_rejection": round(point[tg][1], 4)}
                            for tg in TAGS},
        "class_reweighting": {
            "t25_optimal_for_structural_weight_below": w_hi,
            "t25_optimal_for_structural_weight_above": w_lo,
            "benchmark_structural_fraction": round(frac, 4),
            "note": ("The weight at which 2.5 stops being the pooled argmax, "
                     "compared against the benchmark's own structural fraction. "
                     "Both are computed here; neither is typed.")},
        "paired_bootstrap": {
            "n_draws": N_BOOT, "seed": SEED, "unit": "protein system",
            "recovery": {tg: {"mean": round(float(np.nanmean(rec[tg])), 4),
                              "ci95": ci(rec[tg])} for tg in TAGS},
            "correct_rejection": {tg: {"mean": round(float(np.nanmean(rej[tg])), 4),
                                       "ci95": ci(rej[tg])} for tg in TAGS},
            "p_recovery_fell_t10_to_tSAP": round(p_rec_fell, 4),
            "p_rejection_rose_t10_to_tSAP": round(p_rej_rose, 4)},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canonical", type=Path, default=CANON)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    res = build(pd.read_csv(args.canonical, low_memory=False))

    if args.check:
        if not OUT.exists():
            raise SystemExit("FAIL: %s absent; run without --check first" % OUT.name)
        old = json.loads(OUT.read_text())
        for key in ("population", "point_estimates", "class_reweighting"):
            if old.get(key) != res.get(key):
                raise SystemExit("FAIL: %s changed against the committed record" % key)
        # the bootstrap is seeded, so its summaries must reproduce exactly too
        if old["paired_bootstrap"] != res["paired_bootstrap"]:
            raise SystemExit("FAIL: paired bootstrap changed against the record")
        print("PASS: operating-point diagnostics reproduce the committed record")
        return 0

    p = res["population"]
    print("graded %d = %d structural + %d silent; structural fraction %.4f"
          % (p["n_graded"], p["n_structural"], p["n_silent"],
             p["structural_fraction"]))
    print("\nthreshold        recovery   correct rejection")
    for tg in TAGS:
        e = res["point_estimates"][tg]
        print("  %-13s %.3f      %.3f"
              % (e["threshold_kcal_mol"], e["recovery"], e["correct_rejection"]))
    cw = res["class_reweighting"]
    print("\n2.5 remains the pooled argmax for structural weight <= %s"
          % cw["t25_optimal_for_structural_weight_below"])
    print("  benchmark's own structural fraction: %.4f"
          % cw["benchmark_structural_fraction"])
    pb = res["paired_bootstrap"]
    print("\npaired whole-system bootstrap (%d draws, seed %d, %d systems):"
          % (pb["n_draws"], pb["seed"], p["n_systems"]))
    for tg in (TAGS[0], TAGS[-1]):
        r_, j_ = pb["recovery"][tg], pb["correct_rejection"][tg]
        print("  %-13s recovery %.3f %s   rejection %.3f %s"
              % (KCAL[tg], r_["mean"], r_["ci95"], j_["mean"], j_["ci95"]))
    print("  P(recovery fell) = %.3f   P(rejection rose) = %.3f"
          % (pb["p_recovery_fell_t10_to_tSAP"], pb["p_rejection_rose_t10_to_tSAP"]))

    OUT.write_text(json.dumps(res, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
