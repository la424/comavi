#!/usr/bin/env python3
"""Fetch the structures this example needs. No FoldX licence required.

Downloads the experimental complex 6XI7 from RCSB and splits out the two
single-chain files MAVIS uses for the monomer-fold axis. Chain extraction from
the same crystal (rather than an AlphaFold model) keeps the monomer and complex
axes on identical coordinates, so the fold-vs-complex-fold comparison is not
confounded by a change of structure source.

    python examples/kras_craf/fetch_structures.py --out examples/kras_craf/structures
"""
import argparse
import sys
import urllib.request
from pathlib import Path

RCSB = "https://files.rcsb.org/download"
PDB_ID = "6XI7"
CHAINS = {"A": "fold_kras_model_0.pdb", "B": "fold_raf1_model_0.pdb"}


def fetch(pdb_id, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"{RCSB}/{pdb_id}.pdb"
    print(f"fetching {url}")
    urllib.request.urlretrieve(url, dest)
    if dest.stat().st_size == 0:
        raise SystemExit(f"empty download: {url}")
    print(f"  -> {dest} ({dest.stat().st_size/1024:.0f} kB)")
    return dest


def split_chain(src, chain, dest):
    """Write ATOM/TER records for one chain. HETATM is left out on purpose --
    the pipeline's own preprocessing handles ligand policy per system."""
    kept = 0
    with open(src) as fh, open(dest, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "ANISOU")) and line[21] == chain:
                out.write(line)
                kept += line.startswith("ATOM")
            elif line.startswith("TER") and len(line) > 21 and line[21] == chain:
                out.write(line)
        out.write("END\n")
    if kept == 0:
        raise SystemExit(f"no ATOM records for chain {chain} in {src}")
    print(f"  -> {dest} (chain {chain}, {kept} atoms)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).parent / "structures"))
    args = ap.parse_args()
    out = Path(args.out)

    complex_pdb = fetch(PDB_ID, out / f"{PDB_ID}.pdb")
    for chain, name in CHAINS.items():
        split_chain(complex_pdb, chain, out / name)

    print(f"\nStructures ready in {out}")
    print("Next: python run.py --config examples/kras_craf/systems.yaml \\")
    print(f"        --variants examples/kras_craf/variants.csv --structures {out} \\")
    print("        --out /tmp/mavis_kras --dry-run")
    return 0


if __name__ == "__main__":
    sys.exit(main())
