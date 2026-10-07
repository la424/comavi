#!/usr/bin/env python3
"""Build COMAVI_delta_calibration_points.csv from a curated, cited input.

WHY THIS EXISTS
---------------
This table is the whole basis of the manuscript's calibration section -- the 63
predicted-against-measured comparisons behind Spearman rho, the fitted slope and
the error-structure figure. It is the only place the paper compares predictions
with measurements rather than with other predictors.

Until now it had NO GENERATOR. Nothing in the repository wrote it; six scripts
read it. It entered in a single commit and was hand-maintained thereafter, with
a `source` column carrying three coarse assay-class strings and no per-row
citation. Six of its columns are pure arithmetic on the other two and were
maintained by hand alongside them.

That is the configuration that produces drift, and it hid a real inconsistency.
Four of the five hemoglobin rows come from Kiger/Kwiatkowski, whose own wording
is "9 kcal/tetramer", so they are PER TETRAMER. The fifth, N102T, carried 1.5 --
which is the PER INTERFACE reading of Bonaventura & Riggs, a paper that states
no kcal value at all and reports only a dissociation constant. All five sat in
one table under one axis label and one system label, mixing two conventions.

The convention is now harmonised to per tetramer and enforced here in code:
N102T is DERIVED from Bonaventura's own Kd values rather than typed, so the
convention cannot silently drift back, and the assertion below fails if the
curated input disagrees with that derivation.

INPUT   inputs/calibration_measured_values.csv -- curated, one row per
        comparison, carrying the measured value, the FoldX prediction, and for
        every row a citation and a note on how the measured value was obtained.
OUTPUT  reference_outputs/COMAVI_delta_calibration_points.csv

Usage:
    PYTHONPATH=. python scripts/build_delta_calibration_points.py
    PYTHONPATH=. python scripts/build_delta_calibration_points.py --check
"""
from __future__ import annotations

import argparse
import math
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
SRC = REPO / "inputs" / "calibration_measured_values.csv"
OUT = REPO / "reference_outputs" / "COMAVI_delta_calibration_points.csv"

# The reference threshold the d_meas / d_pred / straddles columns are taken
# relative to. Same 2.5 kcal/mol operating point the manuscript reports.
REF = 2.5

# Gas constant in kcal/(mol K) and the temperature the Kd ratio is read at.
R_GAS = 0.0019872
T_K = 298.15

# Bonaventura & Riggs 1968, JBC 243:980, Hb Kansas. The paper reports
# Kd(tetramer -> dimer) of "about 1.0 to 1.4 x 10^-6 M" for normal
# oxyhemoglobin and "about 200 x 10^-6 M" for Hb Kansas, and states no kcal
# value anywhere. The midpoint of the normal range is used.
KD_NORMAL_UM = (1.0 + 1.4) / 2
KD_KANSAS_UM = 200.0


def n102t_per_tetramer() -> float:
    """N102T assembly penalty, per tetramer, from the source's own Kd values."""
    return round(R_GAS * T_K * math.log(KD_KANSAS_UM / KD_NORMAL_UM), 1)


def build() -> pd.DataFrame:
    cur = pd.read_csv(SRC)
    need = {"variant", "measured_kcal", "foldx_ddg", "axis", "system",
            "source", "in_benchmark", "g73", "citation", "derivation"}
    missing = need - set(cur.columns)
    if missing:
        raise SystemExit("FAIL: curated input lacks columns %s" % sorted(missing))
    if cur.citation.isna().any() or (cur.citation.astype(str).str.strip() == "").any():
        bad = cur.loc[cur.citation.astype(str).str.strip() == "", "variant"].tolist()
        raise SystemExit("FAIL: rows with no citation: %s" % bad[:6])

    # Enforce the per-tetramer convention on the hemoglobin block.
    hb = cur[cur.source == "assembly energetics"]
    want = n102t_per_tetramer()
    got = hb.loc[hb.variant == "N102T", "measured_kcal"]
    if got.empty:
        raise SystemExit("FAIL: N102T absent from the assembly-energetics rows")
    if abs(float(got.iloc[0]) - want) > 0.051:
        raise SystemExit(
            "FAIL: N102T measured value is %.2f but Bonaventura's Kd ratio gives "
            "%.2f kcal/mol per tetramer. The other four hemoglobin rows are per "
            "tetramer (Kiger states '9 kcal/tetramer'), so %.2f would reintroduce "
            "the per-interface/per-tetramer mix this generator exists to remove."
            % (float(got.iloc[0]), want, float(got.iloc[0])))

    df = cur.drop(columns=["citation", "derivation"]).copy()
    df["delta"] = df.foldx_ddg - df.measured_kcal
    df["absdelta"] = df.delta.abs()
    df["straddles"] = (df.measured_kcal >= REF) != (df.foldx_ddg >= REF)
    df["g73"] = df.g73.astype(bool)
    df["d_meas"] = (df.measured_kcal - REF).abs()
    df["d_pred"] = (df.foldx_ddg - REF).abs()
    df["both_near"] = (df.d_meas <= 1.0) & (df.d_pred <= 1.0)
    return df[["variant", "measured_kcal", "foldx_ddg", "axis", "system",
               "source", "in_benchmark", "delta", "absdelta", "straddles",
               "g73", "d_meas", "d_pred", "both_near"]]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="recompute and fail if the committed table differs")
    args = ap.parse_args()
    fresh = build()
    if args.check:
        prev = pd.read_csv(OUT)
        if list(prev.columns) != list(fresh.columns):
            print("FAIL: column set differs")
            return 1
        for col in fresh.columns:
            a, b = prev[col], fresh[col]
            if a.dtype.kind in "fc" or b.dtype.kind in "fc":
                same = pd.Series(a).astype(float).round(4).equals(
                    pd.Series(b).astype(float).round(4))
            else:
                same = a.astype(str).equals(b.astype(str))
            if not same:
                print("FAIL: column %s differs from the committed table" % col)
                return 1
        print("PASS: calibration points reproduce from the curated input "
              "(%d comparisons; N102T at %.1f kcal/mol per tetramer, derived "
              "from Bonaventura's Kd ratio)"
              % (len(fresh), n102t_per_tetramer()))
        return 0
    fresh.to_csv(OUT, index=False)
    print("wrote %s (%d rows)" % (OUT.relative_to(REPO), len(fresh)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
