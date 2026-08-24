from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from build_isds_variant_report import build_outputs  # noqa: E402


class ISDSVariantReportTests(unittest.TestCase):
    def fixture(self) -> pd.DataFrame:
        return pd.DataFrame([
            {
                "system": "CHD_fixture",
                "gene": "GENE1",
                "variant": "A1V",
                "comavi_tier": "Tier 2",
                "ddg_monomer": 1.45,
                "ddg_monomer_confident": True,
                "ddg_fold_partner": 2.90,
                "ddg_binding_partner": 3.50,
                "ddg_partner_confident": True,
                "existing_output": "preserved",
            },
            {
                "system": "CHD_fixture",
                "gene": "GENE2",
                "variant": "G2D",
                "comavi_tier": "Tier 4",
                "ddg_monomer": 0.10,
                "ddg_monomer_confident": True,
                "ddg_fold_partner": 0.20,
                "ddg_binding_partner": 0.30,
                "ddg_partner_confident": True,
                "existing_output": "preserved",
            },
        ])

    def test_historical_csv_is_augmented_and_ranked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "fixture.csv"
            self.fixture().to_csv(source, index=False)
            outputs = build_outputs(source, root / "out", prefix="fixture")
            full = pd.read_csv(outputs["full_csv"])
            ranked = pd.read_csv(outputs["prioritized_csv"])
            summary = json.loads(outputs["summary_json"].read_text(encoding="utf-8"))
            self.assertTrue((full["existing_output"] == "preserved").all())
            self.assertEqual(ranked.iloc[0]["variant"], "A1V")
            self.assertEqual(summary["isds_available_rows"], 2)

    def test_existing_current_fields_are_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "fixture.csv"
            self.fixture().to_csv(source, index=False)
            first = build_outputs(source, root / "first", prefix="fixture")
            second = build_outputs(first["full_csv"], root / "second", prefix="fixture")
            summary = json.loads(second["summary_json"].read_text(encoding="utf-8"))
            self.assertTrue(summary["stored_vs_recomputed_audit"]["all_fields_match"])

    def test_cli_rejects_non_structural_table_without_traceback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "collapsed.csv"

            pd.DataFrame(
                [
                    {
                        "system": "collapsed",
                        "gene": "GENE1",
                        "variant": "A1V",
                    }
                ]
            ).to_csv(source, index=False)

            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        SCRIPTS
                        / "build_isds_variant_report.py"
                    ),
                    str(source),
                    "--out-dir",
                    str(root / "out"),
                ],
                capture_output=True,
                text=True,
            )

            self.assertEqual(
                result.returncode,
                2,
            )

            self.assertIn(
                "No monomer or partner energetic "
                "columns were found",
                result.stderr,
            )

            self.assertNotIn(
                "Traceback",
                result.stderr,
            )


if __name__ == "__main__":
    unittest.main()
