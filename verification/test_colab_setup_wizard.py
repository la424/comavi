from __future__ import annotations

from pathlib import Path
import json
import sys
import tempfile
import unittest

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "notebooks"
if str(NOTEBOOKS) not in sys.path:
    sys.path.insert(0, str(NOTEBOOKS))

import comavi_setup_wizard as wizard

AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


def write_pdb(path: Path, chains: dict[str, tuple[str, list[int]]]) -> None:
    lines = []
    serial = 1
    x = 0.0
    for chain_id, (sequence, numbers) in chains.items():
        if len(sequence) != len(numbers):
            raise ValueError("sequence/number mismatch")
        for amino_acid, residue_number in zip(sequence, numbers):
            lines.append(
                f"ATOM  {serial:5d}  CA  {AA3[amino_acid]:>3s} {chain_id:1s}"
                f"{residue_number:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}"
                f"  1.00{90.0:6.2f}           C  "
            )
            serial += 1
            x += 1.5
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class TestSetupWizard(unittest.TestCase):
    def test_hbb_historical_numbering_is_detected(self) -> None:
        variants = pd.DataFrame(
            [
                {"gene": "HBB", "ref_aa": "E", "position": 6, "alt_aa": "V"},
                {"gene": "HBB", "ref_aa": "E", "position": 6, "alt_aa": "A"},
            ]
        )
        result = wizard.infer_input_offset(
            "HBB",
            variants,
            "MVHLTPEEKSAVTALWGKV",
        )
        self.assertEqual(result.selected_offset, -1)
        self.assertEqual(result.status, "yellow")
        self.assertIn(-2, result.candidate_offsets)
        self.assertEqual(result.matched_variants, 2)

    def test_chain_mapping_and_separate_pipeline_offsets(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "hemoglobin.pdb"
            hbb_reference = "MVHLTPEEKSAVTALWGKV"
            hba_reference = "MVLSPADKTNVKAAWGKV"
            # Mature HBB omits the initiating methionine and starts at raw residue 1.
            write_pdb(
                structure,
                {
                    "A": (hba_reference, list(range(1, len(hba_reference) + 1))),
                    "B": (hbb_reference[1:], list(range(1, len(hbb_reference)))),
                },
            )
            variants = pd.DataFrame(
                [
                    {"gene": "HBB", "ref_aa": "E", "position": 6, "alt_aa": "V"},
                    {"gene": "HBA1", "ref_aa": "D", "position": 7, "alt_aa": "N"},
                ]
            )
            assessment = wizard.assess_structure(
                structure,
                {"hbb": hbb_reference, "hba1": hba_reference},
                variants,
            )
            self.assertEqual(assessment.chain_map["hbb"], "B")
            self.assertEqual(assessment.chain_map["hba1"], "A")
            # HBB historical E6 maps directly to mature-chain residue 6.
            self.assertEqual(assessment.pipeline_offsets["hbb"], 0)
            self.assertEqual(assessment.pipeline_offsets["hba1"], 0)
            self.assertTrue(all(row["covered"] for row in assessment.variant_rows))
            self.assertIn(assessment.status, {"green", "yellow"})

    def test_wrong_chain_override_fails_variant_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "two_chains.pdb"
            write_pdb(
                structure,
                {
                    "A": ("ACDEFGHIK", list(range(1, 10))),
                    "B": ("MMMMMMMMM", list(range(1, 10))),
                },
            )
            variants = pd.DataFrame(
                [{"gene": "GENE1", "ref_aa": "D", "position": 3, "alt_aa": "A"}]
            )
            assessment = wizard.assess_structure(
                structure,
                {"gene1": "ACDEFGHIK"},
                variants,
                chain_overrides={"gene1": "B"},
            )
            self.assertEqual(assessment.status, "red")
            self.assertEqual(assessment.assignment_rows[0].variant_coverage, 0.0)

    def test_nonuniform_structure_numbering_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "nonuniform.pdb"
            write_pdb(
                structure,
                {"A": ("ACDEFG", [1, 2, 3, 10, 11, 12])},
            )
            variants = pd.DataFrame(
                [
                    {"gene": "G", "ref_aa": "C", "position": 2, "alt_aa": "A"},
                    {"gene": "G", "ref_aa": "F", "position": 5, "alt_aa": "A"},
                ]
            )
            assessment = wizard.assess_structure(
                structure,
                {"g": "ACDEFG"},
                variants,
            )
            self.assertEqual(assessment.status, "red")
            self.assertLess(assessment.assignment_rows[0].numbering_uniformity, 0.95)

    def test_preflight_traffic_light(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            monomer = root / "mono.pdb"
            complex_path = root / "complex.pdb"
            sequence = "ACDEFGHIK"
            write_pdb(monomer, {"A": (sequence, list(range(1, 10)))})
            write_pdb(complex_path, {"B": (sequence, list(range(1, 10)))})
            variants = pd.DataFrame(
                [{"gene": "G", "ref_aa": "D", "position": 3, "alt_aa": "A"}]
            )
            mono = wizard.assess_monomer_set(
                {"g": monomer},
                {"g": sequence},
                variants,
            )
            complex_assessment = wizard.assess_structure(
                complex_path,
                {"g": sequence},
                variants,
            )
            preflight = wizard.build_preflight_table(variants, mono, complex_assessment)
            self.assertEqual(preflight.iloc[0]["overall"], "green")
            html = wizard.render_traffic_light_html(preflight)
            self.assertIn("READY", html)

    def test_bundle_roundtrip_and_merge(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "system.pdb"
            write_pdb(structure, {"A": ("ACDE", [1, 2, 3, 4])})
            config = root / "config.yaml"
            config.write_text(
                yaml.safe_dump(
                    {
                        "systems": {
                            "g_system": {
                                "structure_type": "PDB",
                                "pdb_file": "system.pdb",
                                "genes": {
                                    "g": {
                                        "chain": "A",
                                        "monomer_pdb": "system.pdb",
                                        "monomer_offset": 0,
                                    }
                                },
                            }
                        }
                    },
                    sort_keys=False,
                ),
                encoding="utf-8",
            )
            preflight = pd.DataFrame(
                [{"gene": "g", "variant": "D3A", "overall": "green"}]
            )
            bundle = wizard.create_system_setup_bundle(
                root / "bundle.zip",
                config_path=config,
                structure_paths=[structure],
                reference_sequences={"g": "ACDE"},
                preflight=preflight,
                setup_report={"status": "green"},
            )
            manifest = wizard.extract_and_verify_setup_bundle(bundle, root / "extracted")
            self.assertEqual(manifest["setup_wizard_version"], wizard.SETUP_WIZARD_VERSION)
            merged_config, structures, merged_manifest = wizard.merge_setup_bundles(
                [bundle], root / "merged"
            )
            self.assertTrue(merged_config.is_file())
            self.assertTrue((structures / "system.pdb").is_file())
            self.assertEqual(merged_manifest["systems"], ["g_system"])
            self.assertEqual(wizard.genes_from_config(merged_config), ["g"])
            merged_sequences = wizard.read_fasta(root / "merged" / "reference_sequences.fasta")
            self.assertEqual(merged_sequences, {"g": "ACDE"})

    def test_prepared_system_revalidation_rebuilds_offsets_for_current_numbering(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structures = root / "structures"
            structures.mkdir()
            hbb_reference = "MVHLTPEEKSAVTALWGKV"
            hba_reference = "MVLSPADKTNVKAAWGKV"
            hbb_monomer = structures / "hbb_full.pdb"
            hba_monomer = structures / "hba_full.pdb"
            complex_path = structures / "hemoglobin.pdb"
            write_pdb(hbb_monomer, {"A": (hbb_reference, list(range(1, len(hbb_reference) + 1)))})
            write_pdb(hba_monomer, {"A": (hba_reference, list(range(1, len(hba_reference) + 1)))})
            write_pdb(
                complex_path,
                {
                    "A": (hba_reference, list(range(1, len(hba_reference) + 1))),
                    "B": (hbb_reference[1:], list(range(1, len(hbb_reference)))),
                },
            )
            config = root / "config.yaml"
            config.write_text(
                yaml.safe_dump(
                    {
                        "systems": {
                            "hbb_hba1": {
                                "structure_type": "PDB",
                                "pdb_file": complex_path.name,
                                "genes": {
                                    "hbb": {
                                        "chain": "B",
                                        "monomer_pdb": hbb_monomer.name,
                                        "monomer_offset": -1,
                                        "multimer_offset": 0,
                                    },
                                    "hba1": {
                                        "chain": "A",
                                        "monomer_pdb": hba_monomer.name,
                                        "monomer_offset": 0,
                                        "multimer_offset": 0,
                                    },
                                },
                            }
                        }
                    },
                    sort_keys=False,
                ),
                encoding="utf-8",
            )
            # Current batch uses full-reference numbering: HBB Glu7 rather than historical E6.
            variants = pd.DataFrame(
                [{"gene": "HBB", "ref_aa": "E", "position": 7, "alt_aa": "V"}]
            )
            rebuilt, preflight, reports = wizard.revalidate_prepared_systems(
                config,
                structures,
                {"hbb": hbb_reference, "hba1": hba_reference},
                variants,
                output_config_path=root / "config_current.yaml",
            )
            current = yaml.safe_load(rebuilt.read_text(encoding="utf-8"))
            hbb = current["systems"]["hbb_hba1"]["genes"]["hbb"]
            self.assertEqual(hbb["monomer_offset"], 0)
            self.assertEqual(hbb["multimer_offset"], 1)
            self.assertFalse((preflight["overall"] == "red").any())
            self.assertEqual(set(preflight["system"]), {"hbb_hba1"})
            self.assertIn(reports["hbb_hba1"]["status"], {"green", "yellow"})

    def test_config_provenance_is_exact_and_system_scoped(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.yaml"
            config.write_text(
                yaml.safe_dump(
                    {
                        "systems": {
                            "hbb_hba1": {
                                "structure_type": "PDB",
                                "pdb_file": "2HHB.pdb",
                                "genes": {
                                    "hbb": {
                                        "chain": "B",
                                        "monomer_pdb": "hbb.pdb",
                                        "monomer_offset": 0,
                                        "multimer_offset": 1,
                                    },
                                    "hba1": {
                                        "chain": "A",
                                        "monomer_pdb": "hba1.pdb",
                                        "multimer_offset": 1,
                                    },
                                },
                            }
                        }
                    },
                    sort_keys=False,
                ),
                encoding="utf-8",
            )
            single = wizard.build_config_provenance(config)
            self.assertEqual(single["monomer_offsets"], {"hbb": 0, "hba1": 0})
            self.assertEqual(single["multimer_offsets"], {"hbb": 1, "hba1": 1})
            self.assertEqual(single["gene_chain_map"], {"hbb": "B", "hba1": "A"})
            self.assertEqual(single["complex_structure_type"], "PDB")
            self.assertEqual(single["complex_structure"], "2HHB.pdb")

            systems = yaml.safe_load(config.read_text(encoding="utf-8"))["systems"]
            systems["second_system"] = {
                "structure_type": "AF",
                "complex_file": "second.pdb",
                "genes": {
                    "g": {
                        "chain": "C",
                        "monomer_file": "g.pdb",
                        "monomer_offset": -2,
                        "multimer_offset": 3,
                    }
                },
            }
            config.write_text(
                yaml.safe_dump({"systems": systems}, sort_keys=False),
                encoding="utf-8",
            )
            multiple = wizard.build_config_provenance(config)
            self.assertEqual(
                multiple["monomer_offsets"],
                {"hbb_hba1": {"hbb": 0, "hba1": 0}, "second_system": {"g": -2}},
            )
            self.assertEqual(
                multiple["complex_structure"],
                {"hbb_hba1": "2HHB.pdb", "second_system": "second.pdb"},
            )
            self.assertEqual(
                multiple["system_configurations"]["second_system"]["gene_chain_map"],
                {"g": "C"},
            )

    def test_plain_language_cards_do_not_create_probability_claims(self) -> None:
        results = pd.DataFrame(
            [
                {
                    "gene": "g",
                    "variant": "D3A",
                    "isds_v1": 0.72,
                    "isds_dominant_axis": "binding",
                    "comavi_mechanism_t25": "binding_destabilization",
                }
            ]
        )
        frame = wizard.plain_language_summary_frame(results)
        self.assertEqual(frame.iloc[0]["dominant_modeled_effect"], "partner interaction")
        rendered = wizard.render_plain_language_cards(results)
        self.assertIn("not a probability", rendered)
        self.assertIn("interaction assay", rendered)


if __name__ == "__main__":
    unittest.main()
