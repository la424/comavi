#!/usr/bin/env python3
"""Compare a MAVIS run of this example against the frozen benchmark values.

Usage:
    python examples/kras_craf/compare_to_reference.py /tmp/mavis_kras/structural_results.csv

Reports per-variant, per-axis differences against reference_output.csv. FoldX is
stochastic, so exact equality is not expected; the tolerance is set to 3x the
largest per-axis replicate SD recorded in the reference run, which is the scale
of run-to-run variation actually observed rather than an arbitrary epsilon.
"""
import sys
from pathlib import Path

import pandas as pd

REF = Path(__file__).parent / "reference_output.csv"
AXES = [("ddg_monomer", "ddg_monomer_sd"),
        ("ddg_fold_raf1", "ddg_fold_raf1_sd"),
        ("ddg_binding_raf1", "ddg_binding_raf1_sd")]


def main(path):
    ref = pd.read_csv(REF)
    got = pd.read_csv(path, low_memory=False)
    if "variant" not in got.columns:
        got["variant"] = (got.ref_aa.astype(str) + got.position.astype(str)
                          + got.alt_aa.astype(str))

    tol = max(3 * ref[sd].max() for _, sd in AXES if sd in ref.columns)
    print(f"tolerance {tol:.3f} kcal/mol (3x max replicate SD in reference)\n")

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
            print(f"{r.variant:6s} {col:18s} ref {r[col]:8.4f}  got {float(g[col]):8.4f}"
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
