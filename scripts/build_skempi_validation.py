#!/usr/bin/env python3
"""Binding-axis check against SKEMPI 2.0, outside the benchmark cohort.

WHY THIS EXISTS
---------------
The manuscript quotes these statistics (combined Spearman rho, its p value, the
sign-correct count) but until v7.14 they existed only as stdout from the S2 Fig
generator: no committed record held them, so nothing could detect drift. This
script is now the single source; the figure and the prose read the record.

It also records the analysis WITH the barnase Glu73 cluster. All seven Glu73
substitutions are measured as destabilizing but predicted by FoldX as
stabilizing -- a systematic sign inversion at one buried charged site. Since
v7.14 the manuscript leads with all 48 comparisons and reports the fit without
the cluster as a sensitivity view: an exclusion that moves a correlation from
0.34 to 0.55 has to be reported with both numbers, and this one was introduced
together with the analysis, so it is not prespecified and nothing here should
describe it that way.

Run from repo root:
    python scripts/build_skempi_validation.py           # write the record
    python scripts/build_skempi_validation.py --check   # verify, write nothing
"""
import argparse
import json
import math
import pathlib
import sys

import pandas as pd
from scipy.stats import pearsonr, spearmanr

REPO = pathlib.Path(__file__).resolve().parents[1]
BB = REPO / "supplement/skempi/bb_binding_validation.csv"
JT = REPO / "supplement/skempi/jt_binding_validation.csv"
OUT = REPO / "reference_outputs/COMAVI_skempi_validation.json"
GLU73 = r"E[A-Z]?73"  # FoldX mutation codes at barnase Glu73 (chain letter optional)
TOL = 1e-6            # unrounded statistics differ in the last bits across platforms


def stats(df):
    r, pr = pearsonr(df["ddg_meas"], df["ddg_pred_bind"])
    rho, prho = spearmanr(df["ddg_meas"], df["ddg_pred_bind"])
    destab = df[df["ddg_meas"] > 0]
    return {"n": int(len(df)),
            "pearson_r": round(float(r), 4), "pearson_p": float(f"{pr:.4g}"),
            "spearman_rho": round(float(rho), 4), "spearman_p": float(f"{prho:.4g}"),
            "destabilizing_n": int(len(destab)),
            "sign_correct": int((destab["ddg_pred_bind"] > 0).sum())}


def build():
    bb, jt = pd.read_csv(BB), pd.read_csv(JT)
    g73 = bb["foldx_code"].str.match(GLU73)
    cols = ["ddg_meas", "ddg_pred_bind"]
    ex = pd.concat([bb[~g73][cols], jt[cols]])
    inc = pd.concat([bb[cols], jt[cols]])
    cl = bb[g73]
    return {
        "source": ("supplement/skempi/*_binding_validation.csv; SKEMPI 2.0 measured "
                   "binding ddG against FoldX AnalyseComplex predictions"),
        "complexes": {"barnase_barstar_all": stats(bb),
                      "barnase_barstar_excluding_glu73": stats(bb[~g73]),
                      "tem1_blip": stats(jt)},
        "combined_excluding_glu73": stats(ex),
        "combined_including_glu73": stats(inc),
        "glu73_cluster": {
            "n": int(g73.sum()),
            "codes": cl["foldx_code"].tolist(),
            "measured_min": round(float(cl["ddg_meas"].min()), 2),
            "measured_max": round(float(cl["ddg_meas"].max()), 2),
            "predicted_min": round(float(cl["ddg_pred_bind"].min()), 2),
            "predicted_max": round(float(cl["ddg_pred_bind"].max()), 2),
            "all_measured_destabilizing": bool((cl["ddg_meas"] > 0).all()),
            "all_predicted_stabilizing": bool((cl["ddg_pred_bind"] < 0).all()),
            "status": ("single-site systematic sign inversion; included in the "
                       "primary statistics, with the fit without it reported as a "
                       "sensitivity view; introduced with the analysis, not "
                       "prespecified"),
        },
    }


def same(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
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
        print(f"PASS: {OUT.name} reproduces from the committed SKEMPI tables")
        return 0
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    c = out["combined_excluding_glu73"]
    i = out["combined_including_glu73"]
    print(f"wrote {OUT.relative_to(REPO)}: combined rho {c['spearman_rho']:.2f} "
          f"(n={c['n']}), with Glu73 {i['spearman_rho']:.2f} (n={i['n']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
