#!/usr/bin/env python3
"""Pooled predicted-versus-measured calibration statistics.

WHY THIS EXISTS
---------------
COMAVI_delta_calibration_stats.json was committed in one commit together with
the points it summarizes and had no generator, so nothing could tell whether it
still described the points. Eight consumers read it (the calibration figure, the
operating-point builder among them). This script is now its only writer.

It also changes what the headline describes. The committed record fitted 56 of
63 comparisons, leaving out the seven barnase Glu73 substitutions, which are
measured as destabilizing but predicted as stabilizing: a systematic sign
inversion at one buried charged site. That exclusion was introduced together
with the analysis, so it cannot be presented as prespecified. The primary
statistics now use every comparison, and the without-Glu73 fit is reported as a
named sensitivity block. The slope barely moves (0.42 against 0.43); the rank
correlation does (0.48 against 0.63).

The first run of this script reproduced the committed without-Glu73 values to
every stored digit, which validates the method before the primary was switched.

Run from repo root:
    python scripts/build_delta_calibration_stats.py           # write the record
    python scripts/build_delta_calibration_stats.py --check   # verify, write nothing
"""
import argparse
import json
import math
import pathlib
import sys

import pandas as pd
from scipy.stats import linregress, spearmanr

REPO = pathlib.Path(__file__).resolve().parents[1]
POINTS = REPO / "reference_outputs/COMAVI_delta_calibration_points.csv"
OUT = REPO / "reference_outputs/COMAVI_delta_calibration_stats.json"
REFERENCE_T = 2.5
TOL = 1e-6


def sig(x, n):
    return float(f"{x:.{n}g}")


def stats(d, n_all):
    lr = linregress(d["measured_kcal"], d["foldx_ddg"])
    rho = spearmanr(d["measured_kcal"], d["foldx_ddg"])
    delta = d["foldx_ddg"] - d["measured_kcal"]
    rd = spearmanr(d["measured_kcal"], delta)
    st = d["straddles"].fillna(False).astype(bool)
    bn = d["both_near"].fillna(False).astype(bool)
    err = delta.abs()
    return {"n_all": int(n_all), "n_fit": int(len(d)), "n_excl": int(n_all - len(d)),
            "slope": round(float(lr.slope), 4), "intercept": round(float(lr.intercept), 4),
            "rho": round(float(rho.statistic), 4), "p": sig(rho.pvalue, 2),
            "image_of_2p5": round(float(lr.slope * REFERENCE_T + lr.intercept), 2),
            "straddle_n": int(st.sum()), "straddle_borderline": int((st & bn).sum()),
            "straddle_substantive": int((st & ~bn).sum()),
            "med_abs_borderline": round(float(err[st & bn].median()), 2),
            "med_abs_substantive": round(float(err[st & ~bn].median()), 2),
            "rho_delta_vs_measured": round(float(rd.statistic), 4),
            "p_delta_vs_measured": sig(rd.pvalue, 3)}


def build():
    d = pd.read_csv(POINTS)
    g73 = d["g73"].fillna(False).astype(bool)
    out = stats(d, len(d))
    out["population"] = "all comparisons in COMAVI_delta_calibration_points.csv"
    ex = stats(d.loc[~g73], len(d))
    ex["population"] = ("excluding the seven barnase Glu73 substitutions (measured "
                        "destabilizing, predicted stabilizing); introduced with the "
                        "analysis, not prespecified")
    out["sensitivity_excluding_glu73"] = ex
    return out


def same(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, float) or isinstance(b, float):
        return math.isclose(float(a), float(b), rel_tol=TOL, abs_tol=TOL)
    return a == b


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="verify the committed record matches a fresh build; write nothing")
    args = ap.parse_args()
    out = build()
    if args.check:
        if not OUT.exists():
            print(f"FAIL: {OUT.relative_to(REPO)} missing -- run without --check")
            return 1
        if not same(json.loads(OUT.read_text()), out):
            print(f"FAIL: {OUT.name} differs from a fresh build")
            return 1
        print(f"PASS: {OUT.name} reproduces from the committed calibration points")
        return 0
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    ex = out["sensitivity_excluding_glu73"]
    print(f"wrote {OUT.relative_to(REPO)}: all {out['n_fit']} slope {out['slope']:.2f} "
          f"rho {out['rho']:.2f}; without Glu73 {ex['n_fit']} slope {ex['slope']:.2f} "
          f"rho {ex['rho']:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
