#!/usr/bin/env python3
"""Verify the reusable COMAVI setup-wizard module with a synthetic mapping case."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys
import tempfile

import pandas as pd


def write_pdb(path: Path) -> None:
    aa3 = {
        "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
        "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
        "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
        "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    }
    chains = {
        "A": ("MVLSPADKTNVKAAWGKV", list(range(1, 19))),
        "B": ("VHLTPEEKSAVTALWGKV", list(range(1, 19))),
    }
    lines = []
    serial = 1
    coordinate = 0.0
    for chain_id, (sequence, numbers) in chains.items():
        for amino_acid, residue_number in zip(sequence, numbers):
            lines.append(
                f"ATOM  {serial:5d}  CA  {aa3[amino_acid]:>3s} {chain_id:1s}"
                f"{residue_number:4d}    {coordinate:8.3f}{0.0:8.3f}{0.0:8.3f}"
                f"  1.00{90.0:6.2f}           C  "
            )
            serial += 1
            coordinate += 1.5
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--module",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "notebooks" / "comavi_setup_wizard.py",
    )
    args = parser.parse_args()

    module_path = args.module.resolve()
    if not module_path.is_file():
        raise SystemExit(f"SETUP WIZARD CONTRACT: FAIL\nMissing module: {module_path}")
    sys.path.insert(0, str(module_path.parent))
    import comavi_setup_wizard as wizard

    required = {
        "assess_structure",
        "assess_monomer_set",
        "rank_structure_candidates",
        "build_preflight_table",
        "render_traffic_light_html",
        "create_system_setup_bundle",
        "merge_setup_bundles",
        "revalidate_prepared_systems",
        "build_config_provenance",
        "read_fasta",
        "render_plain_language_cards",
    }
    missing = sorted(name for name in required if not hasattr(wizard, name))
    if missing:
        raise SystemExit(f"SETUP WIZARD CONTRACT: FAIL\nMissing public functions: {missing}")
    if wizard.SETUP_WIZARD_VERSION != "COMAVI-Setup-v1":
        raise SystemExit(
            "SETUP WIZARD CONTRACT: FAIL\nUnexpected version: "
            + repr(wizard.SETUP_WIZARD_VERSION)
        )

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        structure = root / "hemoglobin.pdb"
        write_pdb(structure)
        variants = pd.DataFrame(
            [
                {"gene": "HBB", "ref_aa": "E", "position": 6, "alt_aa": "V"},
                {"gene": "HBA1", "ref_aa": "D", "position": 7, "alt_aa": "N"},
            ]
        )
        assessment = wizard.assess_structure(
            structure,
            {
                "hbb": "MVHLTPEEKSAVTALWGKV",
                "hba1": "MVLSPADKTNVKAAWGKV",
            },
            variants,
        )
        if assessment.chain_map != {"hbb": "B", "hba1": "A"}:
            raise SystemExit(
                "SETUP WIZARD CONTRACT: FAIL\nUnexpected chain map: "
                + repr(assessment.chain_map)
            )
        if assessment.pipeline_offsets.get("hbb") != 0:
            raise SystemExit(
                "SETUP WIZARD CONTRACT: FAIL\nHistorical HBB numbering was not reconciled."
            )
        if not all(row["covered"] for row in assessment.variant_rows):
            raise SystemExit(
                "SETUP WIZARD CONTRACT: FAIL\nSynthetic variants were not fully covered."
            )

    print(f"Module: {module_path}")
    print(f"Version: {wizard.SETUP_WIZARD_VERSION}")
    print("Automatic chain assignment: PASS")
    print("Historical-to-reference numbering reconciliation: PASS")
    print("Variant residue coverage: PASS")
    print("Prepared-bundle revalidation and result-card APIs: PASS")
    print("SETUP WIZARD CONTRACT: PASS")


if __name__ == "__main__":
    main()
