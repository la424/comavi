from __future__ import annotations

import unittest

import pandas as pd

from scripts.comavi_v7.isds import (
    ISDS_OUTPUT_COLUMNS,
    add_isds_v1_columns,
    infer_partner_labels,
)


class ISDSOutputContractTests(unittest.TestCase):
    def fixture(self) -> pd.DataFrame:
        return pd.DataFrame([
            {
                "system": "CHD_fixture",
                "variant": "A1V",
                "comavi_tier": "Tier 2",
                "ddg_monomer": 1.45,
                "ddg_monomer_confident": True,
                "ddg_fold_partner": 2.90,
                "ddg_binding_partner": 3.50,
                "ddg_partner_confident": True,
                "existing_output": "preserved",
            }
        ])

    def test_partner_inference(self):
        labels = infer_partner_labels(self.fixture().columns)
        self.assertEqual(labels, ["partner"])

    def test_shared_output_contract(self):
        source = self.fixture()
        out = add_isds_v1_columns(source, ["partner"])
        self.assertEqual(out.loc[0, "existing_output"], "preserved")
        for field in ISDS_OUTPUT_COLUMNS:
            self.assertIn(field, out.columns)
        self.assertTrue(bool(out.loc[0, "isds_available"]))

    def test_batch_invariance_for_external_runner(self):
        source = self.fixture()
        unrelated = source.copy()
        unrelated.loc[0, "variant"] = "G2D"
        unrelated.loc[0, "comavi_tier"] = "Tier 4"
        unrelated.loc[0, "ddg_monomer"] = 0.01
        single = add_isds_v1_columns(source, ["partner"]).loc[0, "isds_v1"]
        batch = add_isds_v1_columns(pd.concat([source, unrelated], ignore_index=True), ["partner"]).loc[0, "isds_v1"]
        self.assertEqual(single, batch)


if __name__ == "__main__":
    unittest.main()
