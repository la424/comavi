#!/usr/bin/env python3
"""
v7.2 canonical-table corrections (three standalone defect fixes).

Operates on reference_outputs/scored_61var_canonical.csv in place (with a
timestamped backup). Each fix is independent, idempotent, and reports its own
before/after so the diff can be audited fix-by-fix.

  FIX 1 — phantom complex-fold ground truth on single-chain systems.
          brca1_brct is a single-chain (tandem-BRCT) system with no partner
          chain, but its `expected_ddg_fold_complex` column was populated as a
          verbatim duplicate of `expected_ddg_monomer` (12/12 rows identical).
          The grader therefore demands agreement on a complex-fold axis that
          cannot physically exist, docking otherwise-correct monomer calls to
          `partial` with missed=fold_complex. Replaced with the explicit
          sentinel `not_applicable`, which the tokenizer maps to not_tested —
          i.e. the axis is excluded from grading rather than silently blank.

  FIX 2 — unpersisted mechanism calls.
          17 rows (brca1_brct 12, cfh_c3b 3, vwf_gpiba 2) carry NO
          `comavi_mechanism_*` call columns. For cfh_c3b and vwf_gpiba the
          *grades* are present and feed the headline, so the shipped table
          could not justify its own numbers for those rows. Re-running the
          shipped classifier reproduces all stored grades exactly (verified
          5 rows x 5 thresholds), proving the calls were computed then
          discarded at the write step. This fix persists them.

  FIX 3 — stale summary grade.
          `mech_consistency_summary` is defined as an alias of
          `mech_consistency_t25` (apply_concordance_v5.py:1594). VHL W117R
          holds `inconsistent` in summary while t25 is NA following the v7.1
          `structurally_uncommitted` correction (ledger Sec. 15). Verification
          and the ledger already carry the corrected headline (MC t2.5 = 0.723,
          graded n = 47); only the released CSV column was left stale. This
          re-derives summary := t25 and recomputes the threshold-stable flag.

NO fix in this script changes a prediction. Fixes 1 and 2 touch only rows that
are currently ungraded or whose grades are reproduced identically; Fix 3 aligns
the CSV with the already-published and already-verified value.

Usage:
    python3 scripts/apply_canonical_corrections_v72.py [--dry-run]
"""
import argparse
import datetime as _dt
import shutil
import sys
from pathlib import Path

import pandas as pd

_THIS_DIR = Path(__file__).resolve().parent
_REPO = _THIS_DIR.parent
sys.path.insert(0, str(_THIS_DIR))

import apply_concordance_v5 as acv5  # noqa: E402

CANONICAL = _REPO / "reference_outputs" / "scored_61var_canonical.csv"

# Systems with no partner chain: the complex-fold and binding axes are
# undefined, not merely unmeasured.
SINGLE_CHAIN_SYSTEMS = {"brca1_brct"}

# Sentinel written into ground-truth axis columns that are undefined for the
# system's topology. Distinguishes "this quantity does not exist here" from
# "nobody measured it" (blank/unknown). Tokenizes to not_tested.
NOT_APPLICABLE = "not_applicable"

MC_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
ALL_TAGS = [tag for tag, _ in acv5.THRESHOLD_SPECS]  # t10 t15 t20 t25 tSAP


def _headline(df, label):
    """Print MC at every threshold INCLUDING the Sapozhnikov operating point."""
    parts = []
    for tag in ALL_TAGS:
        s = df[f"mech_consistency_{tag}"].map(MC_MAP).dropna()
        parts.append(f"{tag} {s.mean():.4f} (n={len(s)})")
    print(f"  {label:9s} " + " | ".join(parts))


def fix1_phantom_axes(df, dry_run=False):
    m = df["system"].isin(SINGLE_CHAIN_SYSTEMS)
    dup = (df.loc[m, "expected_ddg_monomer"].astype(str)
           == df.loc[m, "expected_ddg_fold_complex"].astype(str))
    print(f"FIX 1  single-chain rows={int(m.sum())} "
          f"complex-fold duplicating monomer={int(dup.sum())}")
    already = (df.loc[m, "expected_ddg_fold_complex"].astype(str) == NOT_APPLICABLE)
    if already.all() and m.any():
        print("       already applied — no change")
        return df, 0
    n = int(m.sum())
    if not dry_run:
        df.loc[m, "expected_ddg_fold_complex"] = NOT_APPLICABLE
        df.loc[m, "expected_ddg_binding"] = NOT_APPLICABLE
    print(f"       set expected_ddg_fold_complex and expected_ddg_binding "
          f"-> '{NOT_APPLICABLE}' on {n} rows")
    return df, n


def fix2_persist_calls(df, dry_run=False):
    partners = acv5.discover_partners(df)
    miss = df["comavi_mechanism_t25"].isna()
    print(f"FIX 2  rows missing mechanism calls={int(miss.sum())} "
          f"{dict(df.loc[miss, 'system'].value_counts())}")
    if not miss.any():
        print("       already applied — no change")
        return df, 0

    regraded_mismatch = []
    for idx in df.index[miss]:
        r = df.loc[idx]
        for tag, thr in acv5.THRESHOLD_SPECS:
            call = acv5.classify_mechanism_at(r, partners, thr, tier_col="comavi_tier")
            # Safety: where a grade already exists, the recomputed call must
            # reproduce it, or the write is not a faithful restoration.
            stored = r.get(f"mech_consistency_{tag}")
            if pd.notna(stored):
                g, _, _ = acv5.grade_mechanism_consistency(
                    r, call, r.get("expected_mech_class"),
                    acv5.classify_axis_status(r))
                if str(g) != str(stored):
                    regraded_mismatch.append((r["variant"], tag, g, stored))
            if not dry_run:
                df.at[idx, f"comavi_mechanism_{tag}"] = call
        if not dry_run:
            df.at[idx, "comavi_mechanism"] = df.at[idx, "comavi_mechanism_t25"]

    if regraded_mismatch:
        print("       ABORT — recomputed calls do not reproduce stored grades:")
        for v, t, g, s in regraded_mismatch:
            print(f"         {v} {t}: recomputed={g} stored={s}")
        raise SystemExit(1)
    print(f"       persisted calls for {int(miss.sum())} rows x {len(ALL_TAGS)} "
          f"thresholds; all pre-existing grades reproduced exactly")
    return df, int(miss.sum())


def fix3_stale_summary(df, dry_run=False):
    stale = df["mech_consistency_t25"].astype(str) != df["mech_consistency_summary"].astype(str)
    print(f"FIX 3  rows where summary != t25: {int(stale.sum())} "
          f"{list(df.loc[stale, 'variant'])}")
    if not stale.any():
        print("       already applied — no change")
        return df, 0
    for v in df.loc[stale, "variant"]:
        r = df[df.variant.eq(v)].iloc[0]
        print(f"       {v}: summary {r['mech_consistency_summary']!r} -> "
              f"{r['mech_consistency_t25']!r}")
    if not dry_run:
        df["mech_consistency_summary"] = df["mech_consistency_t25"]
        df["nbhd_mech_consistency_summary"] = df["nbhd_mech_consistency_t25"]

        def _stable(row, prefix=""):
            grades = [row.get(f"{prefix}mech_consistency_{s}") for s in ALL_TAGS]
            return len(set(map(str, grades))) == 1

        df["mech_consistency_threshold_stable"] = df.apply(_stable, axis=1)
        df["nbhd_mech_consistency_threshold_stable"] = df.apply(
            lambda r: _stable(r, prefix="nbhd_"), axis=1)
    return df, int(stale.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--canonical", default=str(CANONICAL))
    args = ap.parse_args()

    path = Path(args.canonical)
    df = pd.read_csv(path, low_memory=False)
    print(f"canonical: {path.name}  rows={len(df)}  cols={len(df.columns)}\n")
    _headline(df, "BEFORE")
    print()

    df, n1 = fix1_phantom_axes(df, args.dry_run)
    df, n2 = fix2_persist_calls(df, args.dry_run)
    df, n3 = fix3_stale_summary(df, args.dry_run)

    print()
    _headline(df, "AFTER")
    s = df["mech_consistency_summary"].map(MC_MAP).dropna()
    print(f"  summary   {s.mean():.4f} (graded n={len(s)})")

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return
    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = path.with_suffix(f".pre_v72_{stamp}.csv")
    shutil.copy(path, backup)
    df.to_csv(path, index=False)
    print(f"\nwrote {path}  (backup {backup.name})")
    print(f"rows changed: fix1={n1} fix2={n2} fix3={n3}")


if __name__ == "__main__":
    main()
