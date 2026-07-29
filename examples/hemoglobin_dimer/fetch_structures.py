#!/usr/bin/env python3
"""Fetch the structures this example needs. No FoldX licence required.

Two sources, because the two axes need different things:

  complex axis -- chains A+B of 2HHB (experimental x-ray, 1.74 A) from RCSB.
                  The alpha1-beta1 interface is measured, not predicted.
  monomer axis -- AlphaFold DB models for HBB (P68871) and HBA1 (P69905),
                  which is what `monomer_pdb` in the config refers to.

The monomer models are NOT interchangeable with the crystal chains. The config's
`monomer_offset: -1` is calibrated to AlphaFold numbering, which retains the
initiator Met that mature globin numbering drops, so CSV position 6 maps to model
residue 7 = Glu. Crystal chain B residue 7 is also a glutamate, so swapping the
sources would score Glu7 instead of Glu6 and produce a plausible-looking but
wrong result. This script keeps them straight.

    python examples/hemoglobin_dimer/fetch_structures.py --out examples/hemoglobin_dimer/structures
"""
import argparse
import sys
import urllib.request
from pathlib import Path

RCSB = "https://files.rcsb.org/download"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"

PDB_ID = "2HHB"
DIMER_CHAINS = ("A", "B")
DIMER_OUT = "2HHB_dimer_AB.pdb"
MONOMERS = {"P68871": ("HBB", "fold_hbb_model_0.pdb"),
            "P69905": ("HBA1", "fold_hba1_model_0.pdb")}
EXPECT_MONOMER = {"fold_hbb_model_0.pdb": (7, "GLU")}   # CSV pos 6 + offset -1


def _get(url, dest):
    urllib.request.urlretrieve(url, dest)
    if dest.stat().st_size == 0:
        raise SystemExit(f"empty download: {url}")
    return dest


def fetch_dimer(out):
    src = out / f"{PDB_ID}.pdb"
    print(f"fetching {RCSB}/{PDB_ID}.pdb")
    _get(f"{RCSB}/{PDB_ID}.pdb", src)
    dest = out / DIMER_OUT
    kept = 0
    with open(src) as fh, open(dest, "w") as o:
        for line in fh:
            if line.startswith(("ATOM", "ANISOU", "HETATM")) and line[21] in DIMER_CHAINS:
                o.write(line)
                kept += line.startswith("ATOM")
            elif line.startswith("TER") and len(line) > 21 and line[21] in DIMER_CHAINS:
                o.write(line)
        o.write("END\n")
    if kept == 0:
        raise SystemExit(f"no ATOM records for chains {DIMER_CHAINS} in {src}")
    print(f"  -> {dest} (chains {'+'.join(DIMER_CHAINS)}, {kept} atoms)")
    # HETATM is kept here on purpose: the config's hetatm_strip_all handles
    # HEM/HOH/PO4 itself, so the pipeline sees the same input the benchmark did.
    return dest


def fetch_monomers(out):
    import json
    made = []
    for acc, (gene, name) in MONOMERS.items():
        meta = json.load(urllib.request.urlopen(f"{AFDB_API}/{acc}", timeout=60))
        if not meta:
            raise SystemExit(f"no AlphaFold DB entry for {acc}")
        url = meta[0].get("pdbUrl")
        if not url:
            raise SystemExit(f"AlphaFold DB has no pdb model for {acc}")
        dest = out / name
        print(f"fetching {gene} ({acc}) {url.rsplit('/', 1)[-1]}")
        _get(url, dest)
        print(f"  -> {dest}")
        made.append(dest)
    return made


def verify(out):
    """Fail loudly if the mutated position does not hold the expected residue."""
    ok = True
    for name, (pos, want) in EXPECT_MONOMER.items():
        p = out / name
        got = None
        for line in open(p):
            if (line.startswith("ATOM") and line[12:16].strip() == "CA"
                    and int(line[22:26]) == pos):
                got = line[17:20].strip()
                break
        status = "ok" if got == want else "WRONG"
        ok &= got == want
        print(f"  {name}: residue {pos} = {got} (expect {want}) {status}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).parent / "structures"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    fetch_dimer(out)
    fetch_monomers(out)
    print("\nverifying residue numbering at the mutated position:")
    if not verify(out):
        raise SystemExit("numbering check FAILED - do not score this run")

    print(f"\nStructures ready in {out}")
    print("Next: python run.py --config examples/hemoglobin_dimer/systems.yaml \\")
    print(f"        --variants examples/hemoglobin_dimer/variants.csv --structures {out} \\")
    print("        --out /tmp/comavi_hb --dry-run")
    return 0


if __name__ == "__main__":
    sys.exit(main())
