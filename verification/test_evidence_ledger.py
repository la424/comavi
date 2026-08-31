from __future__ import annotations

import csv
import json
import shutil
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

from verification.verify_evidence_ledger import (
    CANONICAL_SHA256,
    DEFAULT_CANONICAL,
    DEFAULT_LEDGER,
    DEFAULT_SUMMARY,
    read_csv,
    sha256,
    verify,
)


class EvidenceLedgerVerifierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.ledger_fields, cls.ledger_rows = read_csv(DEFAULT_LEDGER)
        cls.summary = json.loads(DEFAULT_SUMMARY.read_text(encoding="utf-8"))

    def run_case(
        self,
        mutate_rows=None,
        mutate_summary=None,
        mutate_canonical=None,
        mutate_canonical_rows=None,
    ):
        with tempfile.TemporaryDirectory() as temp_directory:
            root = Path(temp_directory)
            ledger = root / "ledger.csv"
            summary = root / "summary.json"
            rows = deepcopy(self.ledger_rows)
            if mutate_rows:
                mutate_rows(rows)
            with ledger.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=self.ledger_fields)
                writer.writeheader()
                writer.writerows(rows)

            summary_data = deepcopy(self.summary)
            if mutate_summary:
                mutate_summary(summary_data)
            summary.write_text(
                json.dumps(summary_data, indent=1) + "\n", encoding="utf-8"
            )

            canonical = DEFAULT_CANONICAL
            if mutate_canonical or mutate_canonical_rows:
                canonical = root / "canonical.csv"
                if mutate_canonical_rows:
                    fields, canonical_rows = read_csv(DEFAULT_CANONICAL)
                    mutate_canonical_rows(canonical_rows)
                    with canonical.open(
                        "w", newline="", encoding="utf-8"
                    ) as handle:
                        writer = csv.DictWriter(handle, fieldnames=fields)
                        writer.writeheader()
                        writer.writerows(canonical_rows)
                else:
                    shutil.copyfile(DEFAULT_CANONICAL, canonical)
                if mutate_canonical:
                    mutate_canonical(canonical)

            return verify(canonical, ledger, summary)

    def test_current_release_passes(self):
        self.assertEqual(
            verify(DEFAULT_CANONICAL, DEFAULT_LEDGER, DEFAULT_SUMMARY), []
        )

    def test_missing_row_fails(self):
        failures = self.run_case(lambda rows: rows.pop())
        self.assertTrue(any("ledger missing" in failure for failure in failures))

    def test_duplicate_row_fails(self):
        failures = self.run_case(
            lambda rows: rows.append(deepcopy(rows[0]))
        )
        self.assertTrue(
            any("duplicate ledger key" in failure for failure in failures)
        )

    def test_token_mismatch_fails(self):
        def mutate(rows):
            rows[0]["expected_token"] = (
                "neutral"
                if rows[0]["expected_token"] != "neutral"
                else "destab"
            )

        failures = self.run_case(mutate)
        self.assertTrue(
            any("expected-token mismatch" in failure for failure in failures)
        )

    def test_blank_metadata_fails(self):
        for field in ("evidence_basis", "evidence_citation"):
            with self.subTest(field=field):
                failures = self.run_case(
                    lambda rows, field=field: rows[0].__setitem__(field, "")
                )
                self.assertTrue(
                    any(f"blank {field}" in failure for failure in failures)
                )

    def test_invalid_vocabularies_fail(self):
        def mutate(rows):
            rows[0]["evidence_type"] = "E0_not_controlled"
            rows[0]["evidence_directness"] = "approximate"

        failures = self.run_case(mutate)
        self.assertTrue(
            any("invalid evidence_type" in failure for failure in failures)
        )
        self.assertTrue(
            any(
                "invalid evidence_directness" in failure
                for failure in failures
            )
        )

    def test_a636p_exact_row_drift_fails(self):
        def mutate(rows):
            for row in rows:
                if (
                    row["system"] == "msh2_msh6"
                    and row["variant"] == "A636P"
                    and row["axis"] == "monomer"
                ):
                    row["evidence_basis"] += " (drift)"
                    return
            raise AssertionError("A636P monomer row absent")

        failures = self.run_case(mutate)
        self.assertTrue(
            any("A636P monomer" in failure for failure in failures)
        )

    def test_summary_drift_fails(self):
        failures = self.run_case(
            mutate_summary=lambda data: data.__setitem__(
                "n_committed_axes", 108
            )
        )
        self.assertTrue(
            any("summary field" in failure for failure in failures)
        )

    def test_canonical_hash_drift_fails(self):
        def mutate(path):
            path.write_bytes(path.read_bytes() + b"\n")

        failures = self.run_case(mutate_canonical=mutate)
        self.assertTrue(
            any("canonical SHA-256 mismatch" in failure for failure in failures)
        )
        self.assertEqual(sha256(DEFAULT_CANONICAL), CANONICAL_SHA256)

    def test_duplicate_canonical_variant_fails(self):
        failures = self.run_case(
            mutate_canonical_rows=lambda rows: rows.append(deepcopy(rows[0]))
        )
        self.assertTrue(
            any(
                "duplicate canonical variant key" in failure
                for failure in failures
            )
        )


if __name__ == "__main__":
    unittest.main()
