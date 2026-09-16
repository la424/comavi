#!/usr/bin/env python3
"""v7.7 canonical repair — propagate the adjudicated mechanism classes into
`axis_signature`.

PROBLEM
-------
`expected_mech_class` is the authoritative ground-truth axis pattern; the grader
consumes it. `axis_signature` is a DERIVED, direction-explicit restatement of the
same information (see derive_axis_signature in apply_concordance_v5.py), and
figure3_axis_competency.py keys its panel b on it, carrying a comment asserting
the two "agree row-for-row".

The v7.7 evidence-ledger correction updated `expected_mech_class` on six rows and
never propagated to `axis_signature`, so that invariant is violated and the figure
plots a stale partition (32/9/6/10 against a canonical that says 29/12/6/10):

    msh2_msh6 C697F    no_structural_effect            -> ppi_destab_mechanism
    kras_craf G12D     no_structural_effect            -> ppi_destab_mechanism
    kras_craf G12V     no_structural_effect            -> ppi_destab_mechanism
    pi3k      E542K    complex_fold_and_binding_destab -> ppi_destab_mechanism
    pi3k      E545K    complex_fold_and_binding_destab -> ppi_destab_mechanism
    mlh1_pms2 L749P    complex_fold_and_binding_destab -> ppi_destab_mechanism

SCOPE — why this is not a wholesale re-derivation
-------------------------------------------------
Re-deriving every row disagrees with the stored column on TWENTY rows, not six.
The other fourteen are pre-existing drift from the v7.2-v7.6 in-place canonical
edits (nine brca1_brct rows whose complex-fold axis was withdrawn by BRCT pooling,
and five cam_cav12/smad4_smad3 rows that re-derive from 'measured negative' to
'never measured'). Those are a separate question with their own history, and nine
of them re-derive to `fold_destab_monomer_only` / `fold_stab_monomer_only` —
labels figure3 has no bucket for, which would silently drop them from panel b and
break its own total-equals-graded assertion.

So this script repairs only rows where the stored `axis_signature` CONTRADICTS the
authoritative `expected_mech_class`. Rows that merely differ from a fresh
re-derivation while remaining consistent with `expected_mech_class` are left alone.

This changes a derived label, not a prediction and not a grade. No ddg_* column and
no mech_consistency_* column is touched; the headline metrics cannot move.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import pathlib
import shutil
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
sys.path.insert(0, str(REPO / "scripts"))
from apply_concordance_v5 import derive_axis_signature  # noqa: E402

# Which mechanism class each axis_signature implies. Mirrors the branch structure
# of derive_axis_signature: binding-positive with a fold axis also positive is
# `mixed_structural`; binding-positive alone is the binding class; fold-only is the
# fold class; nothing positive is silent.
# Only DESTABILIZING signatures and the explicit no-effect signature imply a
# mechanism class unambiguously. Deliberately absent, and therefore never
# asserted against:
#   fold_stab_* / binding_stab_*  — a stabilising expectation is not a
#       destabilising mechanism. brca1_brct R1699Q legitimately carries
#       axis_signature `fold_stab_monomer_and_complex` with expected_mech_class
#       `structurally_silent` (fold-intact curation); recording the direction is
#       the whole reason axis_signature exists, so this is not a contradiction.
#   "uncommitted"                 — asserts nothing about a mechanism.
SIG_IMPLIES = {
    "no_structural_effect": "structurally_silent",
    "monomer_fold_and_binding_destab": "mixed_structural",
    "complex_fold_and_binding_destab": "mixed_structural",
    "binding_destab_fold_intact": "ppi_destab_mechanism",
    "fold_destab_monomer_and_complex": "fold_mechanism",
    "fold_destab_monomer_only": "fold_mechanism",
    "fold_destab_complex_only": "fold_mechanism",
}

# Rows whose expected_mech_class lies outside the four core classes carry a
# curation state (structurally_uncommitted, interface_uncommitted_magnitude, ...)
# that no axis signature encodes. They are skipped rather than "repaired".
CORE_CLASSES = {"structurally_silent", "mixed_structural",
                "ppi_destab_mechanism", "fold_mechanism"}


def contradicting(df: pd.DataFrame) -> pd.Series:
    """Rows whose stored axis_signature contradicts expected_mech_class.

    Restricted to cases where the signature implies a class unambiguously AND
    the expected class is one of the four core classes, so the test cannot fire
    on curation states or stabilising expectations it has no authority over.
    """
    implied = df["axis_signature"].map(SIG_IMPLIES)
    return (implied.notna()
            & df["expected_mech_class"].isin(CORE_CLASSES)
            & (implied != df["expected_mech_class"]))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if any row still contradicts")
    args = ap.parse_args()

    df = pd.read_csv(CANON, low_memory=False)
    bad = contradicting(df)

    if args.check:
        if bad.any():
            print(f"FAIL: {int(bad.sum())} row(s) contradict expected_mech_class")
            print(df.loc[bad, ["system", "variant", "axis_signature",
                               "expected_mech_class"]].to_string(index=False))
            return 1
        print(f"OK: axis_signature consistent with expected_mech_class "
              f"on all {len(df)} rows")
        return 0

    if not bad.any():
        print("nothing to repair; axis_signature already consistent")
        return 0

    fresh = df.apply(derive_axis_signature, axis=1)
    print(f"repairing {int(bad.sum())} contradicting row(s):")
    for i in df.index[bad]:
        print(f"  {df.at[i, 'system']:12} {df.at[i, 'variant']:8} "
              f"{df.at[i, 'axis_signature']:32} -> {fresh.at[i]:32} "
              f"(expected_mech_class {df.at[i, 'expected_mech_class']})")

    # Guard: the repaired value must agree with expected_mech_class.
    for i in df.index[bad]:
        implied = SIG_IMPLIES.get(fresh.at[i])
        if implied != df.at[i, "expected_mech_class"]:
            print(f"ABORT: re-derived {fresh.at[i]!r} implies {implied!r}, "
                  f"not {df.at[i, 'expected_mech_class']!r}")
            return 2

    if args.dry_run:
        print(f"[dry-run] {int(bad.sum())} rows would change; nothing written")
        return 0

    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = CANON.with_name(f"scored_61var_canonical.pre_v77sig_{stamp}.csv")
    shutil.copy(CANON, backup)
    df.loc[bad, "axis_signature"] = fresh[bad]
    df.to_csv(CANON, index=False)
    print(f"wrote {CANON.name}  (backup {backup.name})")
    print(f"rows changed: {int(bad.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
