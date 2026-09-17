#!/usr/bin/env python3
"""v7.8 canonical repair -- restore gene and structure_source on the BRCT rows.

PROBLEM
-------
The 12 `brca1_brct` rows of the released canonical carry NULL `gene` and NULL
`structure_source`. That is 20% of the benchmark and the whole monomer-fold
validation arm, and it is the only place in the table where those two columns
are absent -- every other system has both. The values were never lost from the
repository: `inputs/raw/benchmark_variants_v6.csv` carries `gene = brca1` and
`structure_source = 1JNX.pdb` for all twelve. They were dropped somewhere
between the raw input and the canonical.

The gap surfaced while building the per-variant table that reviewers asked for
(a table cannot report a structure for a fifth of its rows if the canonical
does not hold one), which is the only reason it was ever noticed: no gate reads
these two columns, so nothing failed.

SCOPE
-----
Fill ONLY cells that are (a) currently null AND (b) supplied by the raw input.
Never overwrite a populated cell, never touch a row the raw input does not
cover, and never touch any other column. Rows present in the canonical but
absent from the raw input (5 of 61) are left alone by construction.

No metric reads `gene` or `structure_source`, so the headline numbers must be
bit-identical after this runs; `--check` asserts that separately.

USAGE
    python scripts/apply_brct_provenance_v78.py --check      # exit 1 on drift
    python scripts/apply_brct_provenance_v78.py --dry-run    # report, write nothing
    python scripts/apply_brct_provenance_v78.py              # back up, then write
"""

import argparse
import datetime as _dt
import shutil
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
RAW = REPO / "inputs" / "raw" / "benchmark_variants_v6.csv"
FIELDS = ("gene", "structure_source")

# Columns the headline metrics are computed from. None of them is touched here,
# so they must be identical before and after -- checked, not assumed.
GUARDED = (
    "mech_consistency_t25",
    "expected_mech_class",
    "axis_signature",
    "structural_agreement_n_t25",
    "structural_agreement_d_t25",
)


def variant_key(frame):
    """(system, variant), building `variant` from parts when absent."""
    if "variant" in frame.columns:
        variant = frame["variant"].astype(str)
    else:
        variant = (frame["ref_aa"].astype(str)
                   + frame["position"].astype(str)
                   + frame["alt_aa"].astype(str))
    return list(zip(frame["system"].astype(str), variant))


def pending(canonical, raw):
    """{(system, variant): {field: value}} for cells that are null AND supplied."""
    supply = {}
    raw_keys = variant_key(raw)
    for key, (_, row) in zip(raw_keys, raw.iterrows()):
        supply[key] = dict((f, row[f]) for f in FIELDS if f in raw.columns)

    todo = {}
    for key, (_, row) in zip(variant_key(canonical), canonical.iterrows()):
        if key not in supply:
            continue
        gaps = {}
        for field in FIELDS:
            have = row.get(field)
            give = supply[key].get(field)
            if pd.isna(have) and give is not None and not pd.isna(give):
                gaps[field] = give
        if gaps:
            todo[key] = gaps
    return todo


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if any fillable cell is still null")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change; write nothing")
    args = ap.parse_args()

    canonical = pd.read_csv(CANON, low_memory=False)
    raw = pd.read_csv(RAW)
    todo = pending(canonical, raw)

    if args.check:
        if todo:
            print("FAIL: %d row(s) have a null %s the raw input supplies"
                  % (len(todo), "/".join(FIELDS)))
            for (system, variant), gaps in sorted(todo.items())[:15]:
                print("   %-14s %-8s missing %s" % (system, variant, ", ".join(sorted(gaps))))
            print("\n   remedy: python scripts/apply_brct_provenance_v78.py")
            return 1
        print("PASS: every canonical row covered by the raw input has %s"
              % " and ".join(FIELDS))
        return 0

    if not todo:
        print("nothing to fill; provenance already complete")
        return 0

    systems = sorted(set(system for system, _ in todo))
    print("%d row(s) across %d system(s): %s" % (len(todo), len(systems), ", ".join(systems)))
    for field in FIELDS:
        values = sorted(set(str(gaps[field]) for gaps in todo.values() if field in gaps))
        if values:
            print("   %-18s <- %s" % (field, ", ".join(values)))

    if args.dry_run:
        print("[dry-run] nothing written")
        return 0

    before = canonical[[c for c in GUARDED if c in canonical.columns]].copy()

    keys = variant_key(canonical)
    filled = 0
    for idx, key in zip(canonical.index, keys):
        if key in todo:
            for field, value in todo[key].items():
                canonical.at[idx, field] = value
                filled += 1

    after = canonical[[c for c in GUARDED if c in canonical.columns]]
    if not before.equals(after):
        raise SystemExit("ABORT: a guarded metric column changed; nothing written")

    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = CANON.with_name("scored_61var_canonical.pre_v78prov_%s.csv" % stamp)
    shutil.copy(CANON, backup)
    canonical.to_csv(CANON, index=False)
    print("wrote %s  (backup %s)" % (CANON.name, backup.name))
    print("cells filled: %d across %d row(s); guarded metric columns unchanged" % (filled, len(todo)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
