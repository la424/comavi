#!/usr/bin/env python
"""Characterize the FoldX calling-threshold sweep as a set of operating points.

The sweep across t = 1.0 / 1.5 / 2.0 / 2.5 / tSAP is usually presented as a
robustness check ("the headline holds across thresholds"). It is not one.
Mechanism classes move in OPPOSITE directions across the sweep, so the sweep is
a sensitivity-specificity tradeoff, and each threshold is an operating point
with a distinct competency.

Emits:
  COMAVI_threshold_operating_points.csv       per-threshold sens/spec/balanced
  COMAVI_threshold_class_competency.csv       per-threshold x mechanism class
  COMAVI_threshold_labile_variants.csv        variants whose grade moves
  COMAVI_threshold_operating_points.json      headline stats + severity ladder
"""
import argparse
import json
import pathlib
import sys
import tempfile

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402
import benchmark_population as bp  # noqa: E402

OUTPUTS = ("COMAVI_threshold_operating_points.csv", "COMAVI_threshold_class_competency.csv",
           "COMAVI_threshold_labile_variants.csv", "COMAVI_threshold_operating_points.json")

TAGS = [tag for tag, _ in ac.THRESHOLD_SPECS]
GMAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
from apply_concordance_v5 import (  # noqa: E402
    MECH_STRUCTURAL_CLASSES as _MSC, assert_classes_covered)
# One source of truth: see apply_concordance_v5. The local literal here was
# the fifth independent copy of this list and, like the other four, omitted
# ppi_stab_mechanism -- which silently dropped SMAD4 I500V from the
# structural arm and made sensitivity read 0.696 over 28 rather than 0.672
# over 29. The figure generator's cross-check against the canonical caught it.
STRUCT = set(_MSC)


def measured_image(t, slope, intercept):
    """Invert the measured->FoldX regression: what measured ddG does a FoldX
    threshold t correspond to?"""
    return (t - intercept) / slope


def build(out_dir):
    out = REPO / "reference_outputs"
    df = pd.read_csv(out / "scored_61var_canonical.csv")
    partners = [p for p in ac.discover_partners(df) if f"ddg_{p}_confident" in df.columns]

    # v11: the 56-variant benchmark (benchmark_population), not "every row except
    # the two unobservable BRCT variants". The old selection kept the three
    # candidates with no committed expectation, whose grades are missing, so they
    # were counted as threshold-labile (a missing grade range never equals 0) and
    # n_graded read 59. This file was also never regenerated after v8.2, and the
    # Fig 4 cross-check is what caught it.
    g = bp.benchmark(df).copy()
    for t in TAGS:
        g[f"g_{t}"] = g[f"mech_consistency_{t}"].map(GMAP)
    gcols = [f"g_{t}" for t in TAGS]
    assert g[gcols].notna().all().all(), "every benchmark variant is graded at every threshold"

    # threshold-averaged mechanism consistency, per variant
    g["TA_MC"] = g[gcols].mean(axis=1)
    g["grade_range"] = g[gcols].max(axis=1) - g[gcols].min(axis=1)
    g["threshold_stable"] = g.grade_range == 0

    # committed classes only (drop the two non-committal ground truths)
    gs = g[g.expected_mech_class.isin(STRUCT | {"structurally_silent"})].copy()
    assert len(gs) == len(g), "every benchmark variant carries a committed class"
    gs["truth_struct"] = gs.expected_mech_class.isin(STRUCT)

    # ---- operating points -------------------------------------------------
    cal_b = json.load(open(out / "COMAVI_foldx_scale_calibration.json"))
    cal_p = json.load(open(out / "COMAVI_delta_calibration_stats.json"))

    rows = []
    for tag, spec in ac.THRESHOLD_SPECS:
        tv = spec["fold"] if isinstance(spec, dict) else spec
        sens = gs[gs.truth_struct][f"g_{tag}"].mean()
        spec_ = gs[~gs.truth_struct][f"g_{tag}"].mean()
        rows.append(
            dict(
                threshold=tag,
                foldx_ddg=tv,
                measured_image_benchmark_fit=round(
                    measured_image(tv, cal_b["ols_slope"], cal_b["ols_intercept"]), 2),
                measured_image_pooled_fit=round(
                    measured_image(tv, cal_p["slope"], cal_p["intercept"]), 2),
                n_structural=int(gs.truth_struct.sum()),
                sensitivity=round(float(sens), 3),
                n_silent=int((~gs.truth_struct).sum()),
                specificity=round(float(spec_), 3),
                balanced_MC=round(float((sens + spec_) / 2), 3),
                pooled_MC=round(float(gs[f"g_{tag}"].mean()), 3),
                pooled_MC_full=round(float(g[f"g_{tag}"].mean()), 4),
            )
        )
    ops = pd.DataFrame(rows)

    # ---- per-class competency --------------------------------------------
    cls = (
        g[g.expected_mech_class.isin(STRUCT | {"structurally_silent"})]
        .groupby("expected_mech_class")
        .agg(n=("variant", "size"), **{t: (f"g_{t}", "mean") for t in TAGS},
             threshold_averaged=("TA_MC", "mean"))
        .round(3)
    )

    # ---- threshold-labile variants ---------------------------------------
    lab = gs[~gs.threshold_stable].copy()
    lab["max_abs_ddg"] = lab.apply(lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)
    lab["favours_low"] = lab[["g_t10", "g_t15"]].mean(axis=1) > lab[["g_t25", f"g_{TAGS[-1]}"]].mean(axis=1)
    lab = lab[["variant", "system", "expected_mech_class", "truth_struct",
               "max_abs_ddg", "favours_low"] + gcols].sort_values("max_abs_ddg")

    stats = dict(
        n_graded=int(len(g)),
        n_committed=int(len(gs)),
        class_balance=dict(structural=int(gs.truth_struct.sum()),
                           silent=int((~gs.truth_struct).sum())),
        threshold_averaged_MC=round(float(g.TA_MC.mean()), 4),
        threshold_averaged_balanced_MC=round(
            float((gs[gs.truth_struct].TA_MC.mean() + gs[~gs.truth_struct].TA_MC.mean()) / 2), 3),
        n_threshold_stable=int(g.threshold_stable.sum()),
        n_threshold_labile=int((~g.threshold_stable).sum()),
        pooled_MC_swing=round(float(ops.pooled_MC.max() - ops.pooled_MC.min()), 3),
        balanced_MC_swing=round(float(ops.balanced_MC.max() - ops.balanced_MC.min()), 3),
        labile_split=dict(
            structural_favouring_low=int((lab.truth_struct & lab.favours_low).sum()),
            structural_favouring_high=int((lab.truth_struct & ~lab.favours_low).sum()),
            silent_favouring_low=int((~lab.truth_struct & lab.favours_low).sum()),
            silent_favouring_high=int((~lab.truth_struct & ~lab.favours_low).sum()),
        ),
        median_max_ddg=dict(
            labile_favouring_low=round(float(lab[lab.favours_low].max_abs_ddg.median()), 2),
            labile_favouring_high=round(float(lab[~lab.favours_low].max_abs_ddg.median()), 2),
        ),
        calibration_fits=dict(
            benchmark_n=cal_b["n"], benchmark_slope=cal_b["ols_slope"],
            pooled_n=cal_p["n_fit"], pooled_slope=cal_p["slope"],
        ),
    )

    ops.to_csv(out_dir / "COMAVI_threshold_operating_points.csv", index=False)
    cls.to_csv(out_dir / "COMAVI_threshold_class_competency.csv")
    lab.to_csv(out_dir / "COMAVI_threshold_labile_variants.csv", index=False)
    with open(out_dir / "COMAVI_threshold_operating_points.json", "w") as fh:
        json.dump(stats, fh, indent=1)
    return ops, cls, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="recompute into a temporary directory and fail if any committed output differs")
    args = ap.parse_args()
    out = REPO / "reference_outputs"
    if args.check:
        with tempfile.TemporaryDirectory() as td:
            build(pathlib.Path(td))
            stale = [f for f in OUTPUTS if not (out / f).exists()
                     or (pathlib.Path(td) / f).read_bytes() != (out / f).read_bytes()]
        if stale:
            print("FAIL: operating-point outputs are stale: %s" % stale)
            return 1
        print("PASS: threshold operating points reproduce (%d files)" % len(OUTPUTS))
        return 0
    ops, cls, stats = build(out)
    print(ops.to_string(index=False))
    print()
    print(cls.to_string())
    print()
    print(json.dumps(stats, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
