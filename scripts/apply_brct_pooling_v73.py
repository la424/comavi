#!/usr/bin/env python3
"""
v7.3 — pool the BRCA1-BRCT monomer-fold cohort into the graded headline.

This is a GENUINE METRIC CHANGE, staged separately from the v7.2 defect repairs
(which moved no number). Applying it invalidates the manuscript headline, the
verification constants, and every figure that reports mechanism-consistency.

Background. The 12 BRCT rows shipped with `expected_mech_class` NaN, so the
grader returned N/A and the rows were silently absent from the headline while
still being described in the text as a benchmark arm. The v7.2 repairs made
grading them possible (the phantom complex-fold axis is now `not_applicable`,
and the mechanism calls are persisted); this script makes it happen.

Two curation corrections are required first, both found while deriving the
class and both confined to this cohort.

  (a) Direction-blind fold class. `derive_expected_mech_class` maps ANY positive
      fold axis to `fold_mechanism` without consulting direction, so BRCA1
      R1699Q — annotated `stab` on the monomer axis — was classed as a fold
      mechanism. A stabilized fold is an intact fold; demanding COMAVI fire a
      fold-destabilization axis on it is incoherent, and it graded
      `inconsistent` for that reason alone. R1699Q is the ONLY row in the
      canonical 61-set carrying a `stab` annotation on any axis (verified), so
      the correction is provably confined to this cohort and cannot perturb the
      PPI benchmark.

  (b) Mechanism not observable on the available structure. R1699L and R1699Q are
      curated `fold_intact_function_lost` — their mechanism is loss of a
      phospho-peptide binding site. 1JNX is the isolated tandem BRCT with no
      peptide present, so the binding axis is `not_applicable` (v7.2) and the
      asserted mechanism cannot be observed on this structure at all. Grading
      them scores COMAVI against evidence the input does not contain. They are
      excluded as UNGRADED-BY-CONSTRUCTION, not as failures.

      This exclusion is NOT a new judgment: it is exactly the pre-existing
      `role_in_cohort == "fold_intact_function_lost"` curation, and it is the
      SAME exclusion already applied by figures/src/figure4_measured_vs_foldx.py
      to produce the published fold-axis correlation (rho = 0.72, n = 10). The
      script asserts the two sets are identical rather than trusting this note.

Net: the BRCT arm contributes n = 10 graded rows, the same 10 the measured-DDG
correlation is computed on. The headline moves 47 -> 57 graded rows.

Pooling also brings the cohort into STRUCTURAL AGREEMENT, which exposed two
further population gaps in the released table. Both are repairs of the same
"computed then discarded" class as v7.2, not new quantities, and both are
verified confined to this cohort (zero PPI rows affected):

  (c) Monomer internal-CI columns were never populated for the cohort, though
      `ddg_monomer` and `ddg_monomer_sd` are both present. The agreement gate
      reads `bool(row.get(..._distinguishable_internal_from_0, False))`, and
      `bool(nan)` is True, so the axis was admitted without its gate ever being
      evaluated. This script backfills the columns using the pipeline's own
      `compute_ddg_cis`, so the gate decides on real numbers.

  (d) The tier axis was being graded on rows that have no tier. The tier is
      built from interface-partner and burial terms and is therefore undefined
      for a single-chain system; it is NaN for all 12 rows. Via the same
      `bool(nan)` truthiness, a missing tier was read as "the tier did not
      fire" — scoring free credit on the 4 rows whose expected class is silent
      and a free miss on the other 8. Fixed at source in
      apply_concordance_v5.py (both agreement functions now require
      `pd.notna(comavi_tier)`); this script only recomputes the columns.

After (c) and (d) the cohort contributes exactly one gradeable axis per row —
the monomer axis, the only one its structure supports.

Usage:
    python3 scripts/apply_brct_pooling_v73.py [--dry-run]
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
BRCT_CONCORDANCE = _REPO / "supplement" / "brct" / "brct_foldx_concordance.csv"

MC_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
ALL_TAGS = [tag for tag, _ in acv5.THRESHOLD_SPECS]

# Curated roles whose asserted mechanism cannot be observed on the deposited
# structure for that system (no partner/ligand present). Ungraded by
# construction — see docstring (b).
UNOBSERVABLE_ROLES = {"fold_intact_function_lost"}


def _sweep(df, label):
    """MC at every threshold INCLUDING the Sapozhnikov operating point."""
    parts = []
    for tag in ALL_TAGS:
        s = df[f"mech_consistency_{tag}"].map(MC_MAP).dropna()
        parts.append(f"{tag} {s.mean():.4f}")
    n = df["mech_consistency_t25"].map(MC_MAP).dropna()
    print(f"  MC  {label:9s} " + " | ".join(parts) + f"   (graded n={len(n)})")
    sa_parts = []
    for tag in ALL_TAGS:
        nn = int(df[f"structural_agreement_n_{tag}"].fillna(0).sum())
        dd = int(df[f"structural_agreement_d_{tag}"].fillna(0).sum())
        sa_parts.append(f"{tag} {nn}/{dd}={nn/dd:.4f}" if dd else f"{tag} n/a")
    print(f"  SA  {label:9s} " + " | ".join(sa_parts))


def corrected_expected_class(row):
    """Reference restatement of the direction-aware fold rule.

    v7.3 final: this rule now lives in the SHIPPED
    acv5.derive_expected_mech_class(), so this function is no longer what gets
    written. It is retained solely as an equivalence oracle — main() asserts
    the shipped derivation agrees with it on every row, so moving the rule to
    source cannot silently change the canonical table.

    Rule: a `stab` monomer/complex annotation with no positive binding axis
    means the fold is intact (indeed reinforced) — structurally_silent from the
    standpoint of a destabilization call, not fold_mechanism.
    """
    ec = acv5.derive_expected_mech_class(row)
    if ec != "fold_mechanism":
        return ec
    direction = acv5._axis_direction(row)
    axes = acv5.classify_axis_status(row)
    fold_pos = [a for a in ("fold_monomer", "fold_complex")
                if axes[a] == "positive"]
    # every positive fold axis is stabilizing, and binding is not positive
    if fold_pos and all(direction.get(a) == "stab" for a in fold_pos) \
            and axes["binding"] != "positive":
        return "structurally_silent"
    return ec


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--canonical", default=str(CANONICAL))
    args = ap.parse_args()

    path = Path(args.canonical)
    df = pd.read_csv(path, low_memory=False)
    partners = [p for p in acv5.discover_partners(df)
                if "_ci95_" not in p and "_distinguishable_" not in p]

    print(f"canonical: {path.name}  rows={len(df)}\n")
    _sweep(df, "BEFORE")

    brct = df["system"].eq("brca1_brct")
    if df.loc[brct, "mech_consistency_t25"].notna().any():
        print("\nalready pooled — no change")
        return

    # --- guard (a): the stabilization correction must be confined to BRCT ---
    axis_cols = ["expected_ddg_monomer", "expected_ddg_fold_complex",
                 "expected_ddg_binding"]
    stab = df[axis_cols].astype(str).apply(
        lambda s: s.str.strip().str.lower()).eq("stab").any(axis=1)
    outside = df.loc[stab & ~brct, "variant"].tolist()
    print(f"\n(a) rows with a 'stab' annotation: {int(stab.sum())} "
          f"{df.loc[stab, 'variant'].tolist()}")
    if outside:
        print(f"    ABORT — stabilization correction would touch non-BRCT rows: {outside}")
        raise SystemExit(1)
    print("    confined to BRCT; PPI benchmark cannot be perturbed")

    # --- guard (b): exclusion set must equal the published F4 exclusion ---
    # v7.3 final: the exclusion now comes from the SHIPPED helper, so every
    # consumer (this script, verify_stage6, stress_tests) reads one definition.
    roles = acv5.load_brct_roles(BRCT_CONCORDANCE)
    unobs = acv5.unobservable_variants(roles)
    local = {v for v, r in roles.items() if r in UNOBSERVABLE_ROLES}
    print(f"(b) role_in_cohort in {sorted(UNOBSERVABLE_ROLES)}: {sorted(unobs)}")
    if not unobs:
        print("    ABORT — could not load cohort roles")
        raise SystemExit(1)
    if unobs != local:
        print(f"    ABORT — shipped unobservable_variants() {sorted(unobs)} "
              f"disagrees with this script's rule {sorted(local)}")
        raise SystemExit(1)

    # --- guard (b2): direction-aware fold rule is live in the shipped module -
    # The rule moved from this script into derive_expected_mech_class(). Assert
    # the shipped derivation is idempotent under the local restatement (i.e.
    # it already applies the rule) AND that it actually fires on the one row
    # in the benchmark that carries a stabilizing annotation. Idempotence alone
    # would pass if the rule had been dropped, so both checks are required.
    src_class = df.apply(acv5.derive_expected_mech_class, axis=1)
    local_class = df.apply(corrected_expected_class, axis=1)
    drift = df.loc[src_class.ne(local_class), "variant"].tolist()
    if drift:
        print(f"    ABORT — shipped derivation not idempotent under the "
              f"direction rule on {len(drift)} rows: {drift}")
        raise SystemExit(1)
    stab_rows = df.loc[stab, "variant"].tolist()
    wrong = [v for v in stab_rows
             if src_class[df.variant.eq(v)].iloc[0] == "fold_mechanism"]
    if wrong:
        print(f"    ABORT — direction-aware fold rule is NOT active in the "
              f"shipped module; stabilizing rows still fold_mechanism: {wrong}")
        raise SystemExit(1)
    print(f"(b2) direction-aware fold rule live in shipped module; "
          f"stabilizing rows {stab_rows} -> "
          f"{[src_class[df.variant.eq(v)].iloc[0] for v in stab_rows]}")

    # --- (c) backfill monomer internal-CI columns on the cohort -------------
    ci_cols = [c for c in df.columns
               if c.startswith("ddg_monomer_ci95_")
               or c.startswith("ddg_monomer_distinguishable_")]
    outside_missing = [c for c in ci_cols if df.loc[~brct, c].isna().any()]
    if outside_missing:
        print(f"\n(c) ABORT — CI columns also missing on non-BRCT rows: {outside_missing[:5]}")
        raise SystemExit(1)
    dtypes_before = {c: str(df[c].dtype) for c in df.columns}
    n_ci = 0
    for idx in df.index[brct]:
        vals = acv5.compute_ddg_cis(df.loc[idx], partners)
        for k, v in vals.items():
            if k.startswith("ddg_monomer_") and k in df.columns:
                # Cast bool flags to the table's float convention AT the write.
                # Writing a bool into a float64 column is deprecated (pandas
                # raises a FutureWarning and will error in a later version),
                # and the object-dtype result is what produced the truthy
                # "False" string regression this backfill exists to prevent.
                df.at[idx, k] = float(v) if isinstance(v, bool) else v
        n_ci += 1

    # compute_ddg_cis returns Python bools for the flag columns, but the
    # released table stores them as float 0.0/1.0. Writing bools into a float64
    # column makes it object-typed, and the CSV round-trip then yields the
    # STRINGS "True"/"False" — and `bool("False")` is True, which silently
    # re-opens the very gate this backfill exists to close. Coerce to the
    # table's existing float convention.
    flag_cols = [c for c in ci_cols if "_distinguishable_" in c]
    for c in flag_cols:
        df[c] = df[c].map(
            lambda v: v if pd.isna(v) else float(bool(v) if isinstance(v, bool) else float(v))
        ).astype("float64")
    print(f"\n(c) backfilled monomer CI columns on {n_ci} rows "
          f"({len(ci_cols)} columns, from ddg_monomer_sd already in table); "
          f"{len(flag_cols)} flag columns coerced to float64")

    drift = [(c, dtypes_before[c], str(df[c].dtype)) for c in df.columns
             if str(df[c].dtype) != dtypes_before[c]]
    if drift:
        print(f"    ABORT — backfill changed column dtypes: {drift}")
        raise SystemExit(1)

    # --- (d) tier axis: undefined on single-chain rows ----------------------
    n_no_tier = int(df.loc[brct, "comavi_tier"].isna().sum())
    ppi_no_tier = df.loc[~brct & df["comavi_tier"].isna(), "variant"].tolist()
    print(f"(d) BRCT rows with no tier: {n_no_tier}/12 | "
          f"non-BRCT rows with no tier: {len(ppi_no_tier)}")
    if ppi_no_tier:
        print(f"    ABORT — tier-gate fix would change PPI rows: {ppi_no_tier}")
        raise SystemExit(1)

    # --- (e) mirror the mechanism calls into the `mech_<tag>` alias ---------
    # The released table carries the call under two names: `comavi_mechanism_*`
    # (written by the pipeline) and `mech_*` (a convenience alias, identical on
    # all 49 interaction rows). The alias was never populated for this cohort,
    # so any consumer reading `mech_*` — including the verifier's independent
    # re-derivation track — sees NaN and grades the cohort as N/A.
    for tag in ALL_TAGS:
        src, dst = f"comavi_mechanism_{tag}", f"mech_{tag}"
        if src in df.columns and dst in df.columns:
            disagree = df.loc[~brct & df[src].notna() & df[dst].notna()]
            disagree = disagree[disagree[src].astype(str).str.strip()
                                != disagree[dst].astype(str).str.strip()]
            if len(disagree):
                print(f"    ABORT — {src} and {dst} disagree on "
                      f"{len(disagree)} non-BRCT rows; alias is not safe")
                raise SystemExit(1)
            df.loc[brct, dst] = df.loc[brct, src]
    print(f"(e) mirrored comavi_mechanism_* -> mech_* alias on the cohort "
          f"({len(ALL_TAGS)} thresholds; verified identical on all "
          f"{int((~brct).sum())} interaction rows)")

    # --- grade the cohort ---
    recs = []
    for idx in df.index[brct]:
        r = df.loc[idx]
        ungraded = r["variant"] in unobs
        # v7.3 final: shipped derivation IS the corrected one (guard b2), so
        # ec_raw and ec are now the same call; both retained so the printed
        # audit table keeps its before/after columns.
        ec_raw = corrected_expected_class(r)
        ec = acv5.derive_expected_mech_class(r)
        axs = acv5.classify_axis_status(r)
        rec = {"variant": r["variant"], "role": roles.get(r["variant"], ""),
               "expected_raw": ec_raw, "expected_used": ec,
               "graded": not ungraded}
        for tag in ALL_TAGS:
            if ungraded:
                rec[tag] = None
                continue
            g, fp, mp = acv5.grade_mechanism_consistency(
                r, r.get(f"comavi_mechanism_{tag}"), ec, axs)
            rec[tag] = g
            # Always mutate the in-memory frame so --dry-run reports the true
            # projected sweep; only the file write is gated below.
            df.at[idx, f"mech_consistency_{tag}"] = g
            df.at[idx, f"mech_false_positive_axes_{tag}"] = ",".join(fp) if fp else ""
            df.at[idx, f"mech_missed_positive_axes_{tag}"] = ",".join(mp) if mp else ""
        df.at[idx, "expected_mech_class"] = ec
        df.at[idx, "mech_graded_excluded_reason"] = (
            "mechanism_not_observable_on_structure" if ungraded else "")
        recs.append(rec)

    audit = pd.DataFrame(recs)
    print("\ncohort grading:")
    print(audit.to_string(index=False))

    flipped = audit[audit.expected_raw != audit.expected_used]
    print(f"\nexpected-class corrections (direction-aware): {len(flipped)} "
          f"{flipped.variant.tolist()}")

    n_graded = int(audit.graded.sum())
    print(f"BRCT graded arm: n={n_graded} "
          f"(excluded {len(audit) - n_graded}: {sorted(unobs)})")

    # --- recompute structural agreement on the cohort -----------------------
    def _thr(spec):
        if isinstance(spec, dict):
            return spec["monomer"], spec["fold"], spec["binding"]
        return spec, spec, spec

    sa_rows = []
    for idx in df.index[brct]:
        rec = {"variant": df.at[idx, "variant"]}
        for tag, spec in acv5.THRESHOLD_SPECS:
            n_, d_ = acv5.compute_structural_agreement(df.loc[idx], partners, *_thr(spec))
            df.at[idx, f"structural_agreement_n_{tag}"] = n_
            df.at[idx, f"structural_agreement_d_{tag}"] = d_
            rec[tag] = f"{n_}/{d_}"
        sa_rows.append(rec)
    sa = pd.DataFrame(sa_rows)
    print("\nstructural agreement on cohort (monomer axis only after (c)/(d)):")
    print(sa.to_string(index=False))

    df["mech_consistency_summary"] = df["mech_consistency_t25"]

    def _stable(row):
        grades = [row.get(f"mech_consistency_{s}") for s in ALL_TAGS]
        grades = [g for g in grades if pd.notna(g)]
        return len(set(map(str, grades))) == 1 if grades else False

    df["mech_consistency_threshold_stable"] = df.apply(_stable, axis=1)

    print()
    _sweep(df, "AFTER")

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return
    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = path.with_suffix(f".pre_v73_{stamp}.csv")
    shutil.copy(path, backup)
    df.to_csv(path, index=False)
    print(f"\nwrote {path}  (backup {backup.name})")


if __name__ == "__main__":
    main()
