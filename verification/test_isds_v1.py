from __future__ import annotations

import math
import unittest

import pandas as pd

from scripts.comavi_v7.isds import (
    ISDS_VERSION,
    add_isds_v1_columns,
    calculate_isds_v1,
)


class ISDSV1Tests(unittest.TestCase):
    def row(self, **updates):
        base = {
            "comavi_tier": "Tier 1",
            "ddg_monomer": 0.0,
            "ddg_monomer_confident": True,
            "ddg_fold_partner": float("nan"),
            "ddg_binding_partner": float("nan"),
            "ddg_partner_confident": True,
        }
        base.update(updates)
        return pd.Series(base)

    def test_energy_softsign_anchor_values(self):
        for ratio, expected in [(0, 0), (0.5, 1/3), (1, 0.5), (2, 2/3), (4, 0.8)]:
            row = self.row(ddg_monomer=ratio * 2.9)
            result = calculate_isds_v1(row, ["partner"])
            self.assertAlmostEqual(float(result["isds_energy_component"]), expected, places=7)

    def test_tier_mapping(self):
        expected = {1: 1.0, 2: 2/3, 3: 1/3, 4: 0.0}
        for tier, context in expected.items():
            result = calculate_isds_v1(self.row(comavi_tier=f"Tier {tier}", ddg_monomer=2.9), ["partner"])
            self.assertAlmostEqual(float(result["isds_context_component"]), context, places=7)

    def test_absolute_magnitude_preserves_signed_output(self):
        pos = calculate_isds_v1(self.row(ddg_monomer=2.9), ["partner"])
        neg = calculate_isds_v1(self.row(ddg_monomer=-2.9), ["partner"])
        self.assertEqual(pos["isds_v1"], neg["isds_v1"])
        self.assertEqual(float(neg["isds_dominant_signed_ddg"]), -2.9)

    def test_dominant_axis_and_partner(self):
        result = calculate_isds_v1(self.row(
            ddg_monomer=1.0,
            ddg_fold_partner=2.9,
            ddg_binding_partner=7.0,
        ), ["partner"])
        self.assertEqual(result["isds_dominant_axis"], "binding")
        self.assertEqual(result["isds_dominant_partner"], "partner")
        self.assertAlmostEqual(float(result["isds_energy_ratio_uncapped"]), 2.0, places=7)

    def test_confidence_gate(self):
        result = calculate_isds_v1(self.row(
            ddg_monomer=1.0,
            ddg_fold_partner=29.0,
            ddg_partner_confident=False,
        ), ["partner"])
        self.assertEqual(result["isds_dominant_axis"], "monomer")

    def test_missing_tier_returns_na(self):
        result = calculate_isds_v1(self.row(comavi_tier=None, ddg_monomer=2.9), ["partner"])
        self.assertFalse(bool(result["isds_available"]))
        self.assertTrue(pd.isna(result["isds_v1"]))

    def test_no_evaluable_energy_returns_na(self):
        result = calculate_isds_v1(self.row(ddg_monomer=float("nan")), ["partner"])
        self.assertFalse(bool(result["isds_available"]))
        self.assertTrue(pd.isna(result["isds_v1"]))

    def test_batch_invariance(self):
        a = self.row(ddg_monomer=3.0, comavi_tier="Tier 2")
        b = self.row(ddg_monomer=0.2, comavi_tier="Tier 4")
        single = add_isds_v1_columns(pd.DataFrame([a]), ["partner"]).iloc[0]["isds_v1"]
        batch = add_isds_v1_columns(pd.DataFrame([a, b]), ["partner"]).iloc[0]["isds_v1"]
        self.assertEqual(single, batch)

    def test_monotonicity(self):
        low = calculate_isds_v1(self.row(ddg_monomer=1.0, comavi_tier="Tier 2"), ["partner"])
        high_energy = calculate_isds_v1(self.row(ddg_monomer=3.0, comavi_tier="Tier 2"), ["partner"])
        high_context = calculate_isds_v1(self.row(ddg_monomer=1.0, comavi_tier="Tier 1"), ["partner"])
        self.assertGreater(float(high_energy["isds_v1"]), float(low["isds_v1"]))
        self.assertGreater(float(high_context["isds_v1"]), float(low["isds_v1"]))

    def test_version(self):
        result = calculate_isds_v1(self.row(ddg_monomer=1.0), ["partner"])
        self.assertEqual(result["isds_version"], ISDS_VERSION)


if __name__ == "__main__":
    unittest.main()
