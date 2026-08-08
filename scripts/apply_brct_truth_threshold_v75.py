#!/usr/bin/env python3
"""
v7.5 canonical-table correction — align BRCT fold ground truth with the
threshold the pipeline is graded at.

PROBLEM
-------
Every BRCA1 tandem-BRCT variant carrying a positive measured GdmCl unfolding
free energy (Rowling 2010) was annotated `expected_ddg_monomer = destab` and
`role_in_cohort = monomer_fold_destabilizer`, regardless of magnitude. But the
pipeline is graded at a calling threshold of 2.5 kcal/mol. Two variants sit
below it:

    R1751Q   measured dG_U-F = 1.57   FoldX monomer = 0.73   ClinVar Benign
    L1664P   measured dG_U-F = 1.18   FoldX monomer = 0.53   ClinVar Benign

FoldX called both not-destabilizing, which AGREES with the measurement, and the
grader scored both `inconsistent` for it. The annotation held the pipeline to a
standard its own ground truth does not support: a variant destabilized by
1.2-1.6 kcal/mol is not a destabilizer under a 2.5 kcal/mol criterion.

FIX
---
Re-derive the BRCT monomer-fold expectation from the measured value at the same
threshold the prediction is judged at:

    expected_ddg_monomer := 'destab' if measured_dG_U-F >= 2.5 else 'neutral'

Applied ONLY to brca1_brct rows, which are the only rows in the benchmark with a
directly measured subunit stability to re-derive from. Variants whose
expectation flips to neutral and which have no other structural expectation fall
to `expected_mech_class = structurally_silent`, consistent with how the two
measured neutral controls (M1663K -0.03, P1806A 0.06) are already handled.

R1699L/Q are NOT touched: role_in_cohort `fold_intact_function_lost`, measured
-0.99 and -1.83 (stabilizing), already expected-neutral on the fold axis with
their pathogenicity carried by a documented binding-site mechanism.

This changes GROUND TRUTH, not a prediction. No ddg_* column is modified. The
grades that move do so because the standard was corrected, and every affected
row is reported before/after so the diff can be audited variant-by-variant.

STATUS: TESTED AND REJECTED -- DO NOT RUN
-----------------------------------------
This script is retained as the record of a considered alternative. It was
evaluated as a counterfactual (never written to the canonical table) and the
sign-based convention was retained. Do not re-litigate without reading this.

The premise above is wrong in a specific way. It assumes the 2.5 kcal/mol
calling threshold means the same thing on a measured free energy as it does on
a FoldX DDG. It does not. Across the 15 benchmark variants with a directly
measured comparator (10 BRCT fold + 5 hemoglobin binding), FoldX tracks the
measurement in RANK (Spearman rho = 0.81, p = 2.2e-4) but COMPRESSES it in
SCALE (OLS slope 0.61). A measured 2.5 kcal/mol maps to ~2.16 FoldX units. The
sign-based convention is what the annotation uses precisely because the two
scales share a direction but not a unit.

Applied honestly to ALL measured variants -- not only the two that flatter the
fold arm -- six expectations flip and five grades move:

    fixed:  R1751Q (1.57), L1664P (1.18), W37Y (2.00)
    broken: V1808A (2.40), V1665M (2.22)   <- FoldX 3.57 / 2.98, correct today

    headline MC   0.7193 -> 0.7281   (delta +0.0088, n = 57)
    fold class    0.550 (n=10) -> 0.583 (n=6)
    silent class  0.812 (n=32) -> 0.778 (n=36)

That is +0.009 on the headline in exchange for a 40% smaller fold arm -- the
study's weakest and the one the planned BRCT expansion is designed to enlarge --
and a degraded silent class. Rejected on that trade.

The reporting decision instead: keep the convention, report the fold arm at
0.55 (n = 10), and name its four errors explicitly (V1736A and R755W as genuine
misses; R1751Q and L1664P as sub-threshold measured destabilizations where
FoldX agrees with the physics). The scale-compression result is reported in the
Discussion as the reason the conventions differ. See
docs/COMAVI_section_fold_arm_and_calibration.md.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import shutil
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
BRCT = REPO / "supplement" / "brct" / "brct_foldx_concordance.csv"

SYSTEM = "brca1_brct"
THRESHOLD = 2.5
MC_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}


def _headline(df: pd.DataFrame, tag: str) -> None:
    s = df["mech_consistency_t25"].map(MC_MAP).dropna()
    fold = df[df["expected_mech_class"].eq("fold_mechanism")]
    fs = fold["mech_consistency_t25"].map(MC_MAP).dropna()
    print(f"  [{tag}] MC t2.5 = {s.mean():.4f} (graded n={len(s)}) | "
          f"fold class MC = {fs.mean():.4f} (n={len(fs)})")


def retruth(df: pd.DataFrame, dry: bool) -> tuple[pd.DataFrame, int]:
    meas = (pd.read_csv(BRCT)
              .set_index("variant")["measured_ddG_UF_kcal_mol"])

    m = df["system"].eq(SYSTEM)
    if not m.any():
        print("  no brca1_brct rows; nothing to do")
        return df, 0

    changed = 0
    for i in df[m].index:
        v = df.at[i, "variant"]
        if v not in meas.index or pd.isna(meas[v]):
            continue
        # R1699L/Q carry a documented fold-intact/function-lost mechanism.
        if str(df.at[i, "curated_mechanism"]) == "binding_site_loss_fold_intact":
            continue
        want = "destab" if meas[v] >= THRESHOLD else "neutral"
        have = str(df.at[i, "expected_ddg_monomer"])
        if have == want:
            continue
        print(f"    {v:8s} measured {meas[v]:5.2f}  "
              f"expected_ddg_monomer {have!r} -> {want!r}")
        changed += 1
        if not dry:
            df.at[i, "expected_ddg_monomer"] = want
            if want == "neutral":
                df.at[i, "curated_mechanism"] = "no_structural_effect"
                df.at[i, "expected_mech_class"] = "structurally_silent"
    return df, changed


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    df = pd.read_csv(CANON, low_memory=False)
    _headline(df, "BEFORE")
    print(f"\n  re-deriving BRCT fold expectation at {THRESHOLD} kcal/mol:")
    df, n = retruth(df, args.dry_run)
    print()

    if args.dry_run:
        print(f"[dry-run] {n} rows would change; nothing written")
        return

    _headline(df, "AFTER ")
    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = CANON.with_suffix(f".pre_v75_{stamp}.csv")
    shutil.copy(CANON, backup)
    df.to_csv(CANON, index=False)
    print(f"\nwrote {CANON.name}  (backup {backup.name})")
    print(f"rows changed: {n}")
    print("\nNOTE: grades are NOT re-derived here. Re-run the concordance "
          "pipeline to regrade against the corrected ground truth.")


if __name__ == "__main__":
    main()
