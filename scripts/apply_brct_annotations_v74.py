#!/usr/bin/env python3
"""
v7.4 — backfill the structural-provenance annotations the BRCA1-BRCT merge dropped.

This is a DEFECT REPAIR, not a metric change. It moves no headline number
(verified by assertion below): mechanism-consistency and structural agreement
gate on the per-axis `*_confident` / `*_distinguishable_*` flags, which the BRCT
rows already carry. What it repairs is the FOUR-WAY CONCORDANCE FRAMEWORK, which
gates on `ddg_confidence_derived` and was silently dropping all 12 cohort rows.

Background. `site_plddt_status`, `monomer_plddt`, `best_plddt` and
`structure_evaluable` are curated structural-provenance annotations carried on
the PPI intermediate table (`inputs/intermediate/comavi_v7_results_with_nbhd.csv`).
The BRCA1-BRCT cohort was scored in a SEPARATE supplement run
(`supplement/brct/`) and merged in without them, so all 12 rows carry NaN.

Consequence. `derive_ddg_confidence()` short-circuits to "high" on
`site_plddt_status == "crystal"`; failing that it reads `monomer_plddt`, which
is also NaN, and falls through to "low". `axis_evaluable()` then admits the ddg
axis only for confidence in ("high", "medium"), so every BRCT row was excluded
from the concordance denominator — 0/12 rows carry a concordance value.

This is NOT a pipeline defect. No code in `scripts/comavi_v7/` writes
`site_plddt_status`; the pipeline handles crystal structures correctly via
`MultimerSpec.plddt_gate=False`, which makes `get_plddt()` return the sentinel
100.0 everywhere (structure_loading.py). The gap is confined to the benchmark
table's merge of a separately-scored cohort.

The correct values are not a judgment call. PDB 1JNX is an X-ray crystal
structure of the BRCA1 tandem BRCT domain, so it takes exactly the annotation
every other experimental structure in the benchmark takes (2HHB, 6XI7, 1JM7):
`site_plddt_status = "crystal"`, which is the same bypass the pipeline applies
to crystal/NMR inputs. `structure_evaluable = True` follows — the fold axis is
scored on this structure throughout, including for the published rho = 0.72
correlation, so treating it as unevaluable is self-contradictory.

`monomer_plddt` / `best_plddt` are deliberately left NaN: they are pLDDT values
and a crystal structure has none. The "crystal" short-circuit is what makes them
unnecessary, which is precisely how the other experimental systems behave.

SECOND DEFECT, found while verifying this one. Recomputing the concordance
columns also changes 5 NON-cohort rows: cfh R78G / R53H / I62V and vwf R1334Q /
A1381T. Their stored concordance was computed during the v7 expansion BEFORE
their AlphaMissense and clinical annotations were merged, so it was frozen at a
2-axis denominator (e.g. R78G "2/2") while the rows now support 4 axes ("3/4";
the two VWF rows support 3, lacking an AlphaMissense score — the mature protein
exceeds the AM release length limit). Every one of the 5 is `axis_evaluable`
True on struct/ddg/franklin today. Recomputing from the shipped
`compute_concordance()` brings them onto the same denominator as the other 44
PPI rows. Also repaired: `concordance_pathonly_*` was populated on only 5 of 61
rows for the same reason, and is now computed benchmark-wide.

Neither repair moves mechanism-consistency or structural agreement.

Run:  python3 scripts/apply_brct_annotations_v74.py [--dry-run]
"""
import argparse
import datetime as _dt
import pathlib
import shutil
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
SYSTEM = "brca1_brct"

# The annotations to backfill, and the value each takes for an X-ray structure.
BACKFILL = {
    "site_plddt_status": "crystal",
    "structure_evaluable": True,
}

# Non-cohort rows whose stored concordance predates the merge of their external
# annotations (see "SECOND DEFECT" above). The repair refreshes these; any row
# outside this set changing is a regression, not a repair.
EXPECTED_STALE_CONCORDANCE = {
    ("cfh", "R78G"), ("cfh", "R53H"), ("cfh", "I62V"),
    ("vwf", "R1334Q"), ("vwf", "A1381T"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--canonical", default=str(CANON))
    args = ap.parse_args()

    sys.path.insert(0, str(REPO / "scripts"))
    import apply_concordance_v5 as ac

    path = pathlib.Path(args.canonical)
    df = pd.read_csv(path, low_memory=False)
    partners = [p for p in ac.discover_partners(df) if p]
    mask = df["system"].eq(SYSTEM)
    print(f"{SYSTEM}: {int(mask.sum())} rows")

    for col in BACKFILL:
        n_na = int(df.loc[mask, col].isna().sum())
        print(f"  {col:24s} NaN on {n_na}/{int(mask.sum())} cohort rows")

    # --- guard: this repair must be confined to the cohort ---------------
    for col in BACKFILL:
        off = df.loc[~mask, col].isna().sum()
        assert off == 0, f"{col} also NaN on {off} non-cohort rows — repair is not confined"

    # --- baseline metrics BEFORE, to prove nothing moves -----------------
    MCM = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}

    def headline(frame):
        n = d = 0
        for _, r in frame.iterrows():
            for _k, (a, b) in ac.structural_agreement_by_axis(
                    r, partners, 2.5, 2.5, 2.5).items():
                n += a
                d += b
        mc = frame["mech_consistency_t25"].map(MCM)
        return round(mc.mean(), 4), int(mc.notna().sum()), n, d

    before = headline(df)
    conc_before = int(df["concordance_pathonly_full"].notna().sum())

    # --- apply ------------------------------------------------------------
    dtypes_before = {c: str(df[c].dtype) for c in df.columns}
    conc_before_vals = {c: df[c].copy() for c in df.columns if "concordance" in c}
    for col, val in BACKFILL.items():
        # These columns are all-NaN float64 on read when the cohort is the only
        # gap; cast explicitly so assigning a str/bool is not a silent upcast.
        df[col] = df[col].astype(object)
        df.loc[mask, col] = val

    df["ddg_confidence_derived"] = df.apply(
        lambda r: ac.derive_ddg_confidence(r, partners), axis=1)

    cohort_conf = set(df.loc[mask, "ddg_confidence_derived"])
    print(f"  cohort ddg_confidence_derived after repair: {cohort_conf}")
    assert cohort_conf == {"high"}, f"expected all-high on a crystal structure, got {cohort_conf}"

    # --- recompute the concordance columns the repair unblocks ------------
    # Mirrors the Stage-5 call pattern in apply_concordance_v5.py (both
    # pipelines x both modes) so the columns are produced by the shipped
    # implementation rather than a parallel one here.
    specs = [
        ("comavi_tier", "", True),
        ("comavi_tier", "", False),
        ("nbhd_tier", "nbhd_", True),
        ("nbhd_tier", "nbhd_", False),
    ]
    for tier_col, prefix, include_external in specs:
        src = df.apply(
            lambda r, _tc=tier_col, _p=prefix, _ie=include_external:
                ac.compute_concordance(r, r[_tc], prefix=_p, include_external=_ie),
            axis=1, result_type="expand")
        for c in src.columns:
            df[c] = src[c].values

    # Stage 5b: the structural-signal / external-consensus split, same pattern.
    for tier_col, prefix in (("comavi_tier", ""), ("nbhd_tier", "nbhd_")):
        src = df.apply(
            lambda r, _tc=tier_col, _p=prefix:
                ac.compute_signal_consensus_split(r, r[_tc], prefix=_p),
            axis=1, result_type="expand")
        for c in src.columns:
            if prefix == "nbhd_" and "external_consensus" in c:
                continue  # pipeline-independent; taken from p1 only
            df[c] = src[c].values

    n_conc = int(df.loc[mask, "concordance_pathonly_full"].notna().sum())
    print(f"  cohort rows now carrying a concordance value: "
          f"{conc_before} -> {int(df['concordance_pathonly_full'].notna().sum())} overall, "
          f"{n_conc}/{int(mask.sum())} in cohort")

    # --- guard: headline must NOT move ------------------------------------
    after = headline(df)
    assert before == after, f"headline moved: {before} -> {after} (this is a defect repair, not a metric change)"
    print(f"  headline unchanged: MC {after[0]} (n={after[1]}), SA {after[2]}/{after[3]}")

    # --- report: non-cohort rows whose stale concordance is refreshed -------
    # Expected and documented in the module docstring (second defect): rows
    # whose stored concordance predates the merge of their external
    # annotations. Listed explicitly so the change is never silent.
    changed = []
    for c, old in conc_before_vals.items():
        if c.endswith(("_n", "_denom")):
            continue
        new = df[c]
        both = old.notna() & ~mask
        for i in df.index[both & (old.astype(str) != new.astype(str))]:
            changed.append((df.at[i, "gene"], df.at[i, "variant"], c,
                            old.at[i], new.at[i]))
    if changed:
        print(f"  refreshed stale concordance on {len({(g, v) for g, v, *_ in changed})} "
              f"non-cohort variants:")
        for g, v, c, o, n in sorted(changed):
            print(f"    {g:6s} {v:8s} {c:34s} {o} -> {n}")
    assert {(g, v) for g, v, *_ in changed} <= EXPECTED_STALE_CONCORDANCE, (
        f"unexpected non-cohort concordance change: "
        f"{ {(g, v) for g, v, *_ in changed} - EXPECTED_STALE_CONCORDANCE }")

    # --- guard: no dtype drift on untouched columns ------------------------
    # The concordance count columns legitimately go float64 -> int64: they were
    # float only because the stale rows left NaNs, which the recompute fills.
    drift = [(c, dtypes_before[c], str(df[c].dtype)) for c in df.columns
             if c not in BACKFILL and c != "ddg_confidence_derived"
             and dtypes_before[c] != str(df[c].dtype)
             and not (c.endswith(("_n", "_denom"))
                      and any(k in c for k in ("concordance", "structural_signal",
                                               "external_consensus"))
                      and dtypes_before[c] == "float64"
                      and str(df[c].dtype) == "int64")]
    assert not drift, f"unexpected dtype drift: {drift}"

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return

    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = path.with_suffix(f".pre_v74_{stamp}.csv")
    shutil.copy(path, backup)
    df.to_csv(path, index=False)
    print(f"\nwrote {path}  (backup {backup.name})")


if __name__ == "__main__":
    main()
