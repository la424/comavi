#!/usr/bin/env python3
"""Compare a MAVIS run of this example against the frozen benchmark values.

Usage:
    python examples/hemoglobin_dimer/compare_to_reference.py /tmp/mavis_hb/structural_results.csv

FoldX is stochastic, so exact equality is not expected. The tolerance is 3x the
largest per-axis replicate SD recorded in the reference run -- the scale of
run-to-run variation actually observed, rather than an arbitrary epsilon. There
is a floor of 0.15 kcal/mol because one axis here has SD 0.0 (a single-run
value), and a zero tolerance would flag ordinary noise as disagreement.

Mechanism strings are compared exactly: a mechanism flip is a real difference
even when the underlying dG moved only a little.
"""
import sys
from pathlib import Path

import pandas as pd

REF = Path(__file__).parent / "reference_output.csv"
AXES = [("ddg_monomer", "ddg_monomer_sd"),
        ("ddg_fold_hba1", "ddg_fold_hba1_sd"),
        ("ddg_binding_hba1", "ddg_binding_hba1_sd")]
TOL_FLOOR = 0.15


def main(path):
    ref = pd.read_csv(REF)
    got = pd.read_csv(path, low_memory=False)
    if "variant" not in got.columns:
        got["variant"] = (got.ref_aa.astype(str) + got.position.astype(str)
                          + got.alt_aa.astype(str))

    sds = [ref[sd].max() for _, sd in AXES if sd in ref.columns]
    tol = max(TOL_FLOOR, 3 * max(sds)) if sds else TOL_FLOOR
    print(f"tolerance {tol:.3f} kcal/mol (max of {TOL_FLOOR} floor, 3x max replicate SD)\n")

    fails = 0
    for _, r in ref.iterrows():
        g = got[got.variant == r.variant]
        if g.empty:
            print(f"{r.variant:6s} MISSING from run output")
            fails += 1
            continue
        g = g.iloc[0]
        for col, _sd in AXES:
            if col not in got.columns:
                print(f"{r.variant:6s} {col:18s} column absent from run output")
                fails += 1
                continue
            d = float(g[col]) - float(r[col])
            ok = abs(d) <= tol
            fails += not ok
            print(f"{r.variant:6s} {col:18s} ref {float(r[col]):8.4f}  got {float(g[col]):8.4f}"
                  f"  d {d:+7.4f}  {'ok' if ok else 'DIFFERS'}")
        if "mavis_mechanism" in got.columns:
            same = str(g["mavis_mechanism"]) == str(r["mavis_mechanism"])
            fails += not same
            print(f"{r.variant:6s} {'mechanism':18s} ref {r['mavis_mechanism']!r}"
                  f"  got {g['mavis_mechanism']!r}  {'ok' if same else 'DIFFERS'}")
        print()

    print("PASS - reproduces the benchmark values" if not fails
          else f"{fails} difference(s) outside tolerance")
    return 1 if fails else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
