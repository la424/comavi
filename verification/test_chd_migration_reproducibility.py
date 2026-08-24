from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]

SCRIPT = (
    REPO
    / "verification"
    / "verify_chd_migration_reproducibility.py"
)

CSV_FILES = (
    "chd_concordance_with_isds_v1.csv",
    "chd_concordance_prioritized.csv",
)

JSON_FILE = (
    "chd_concordance_isds_summary.json"
)

MARKDOWN_FILE = (
    "chd_concordance_isds_report.md"
)


class TestCHDMigrationReproducibility(
    unittest.TestCase
):
    def write_bundle(
        self,
        root: Path,
        *,
        numeric: str,
        text: str,
    ) -> None:
        root.mkdir(
            parents=True,
            exist_ok=True,
        )

        csv_text = (
            "system,variant,x,label\n"
            f"system_a,V1,{numeric},{text}\n"
        )

        for name in CSV_FILES:
            (root / name).write_text(
                csv_text,
                encoding="utf-8",
            )

        (root / JSON_FILE).write_text(
            json.dumps(
                {
                    "rows": 1,
                    "isds_version": "ISDS-v1",
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

        (root / MARKDOWN_FILE).write_text(
            "# Report\n",
            encoding="utf-8",
        )

    def run_comparison(
        self,
        expected: Path,
        observed: Path,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--expected-dir",
                str(expected),
                "--observed-dir",
                str(observed),
            ],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_equivalent_numeric_serialization_passes(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = root / "expected"
            observed = root / "observed"

            self.write_bundle(
                expected,
                numeric="0.30000000000000004",
                text="same",
            )

            self.write_bundle(
                observed,
                numeric="0.3",
                text="same",
            )

            result = self.run_comparison(
                expected,
                observed,
            )

            self.assertEqual(
                result.returncode,
                0,
                result.stdout + result.stderr,
            )

            self.assertIn(
                "byte_identical=False",
                result.stdout,
            )

            self.assertIn(
                "CHD MIGRATION REPRODUCIBILITY: PASS",
                result.stdout,
            )

    def test_real_numeric_difference_fails(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = root / "expected"
            observed = root / "observed"

            self.write_bundle(
                expected,
                numeric="0.30",
                text="same",
            )

            self.write_bundle(
                observed,
                numeric="0.31",
                text="same",
            )

            result = self.run_comparison(
                expected,
                observed,
            )

            self.assertNotEqual(
                result.returncode,
                0,
            )

            combined = (
                result.stdout
                + result.stderr
            )

            self.assertIn(
                "column=x",
                combined,
            )

            self.assertIn(
                "CHD MIGRATION REPRODUCIBILITY: FAIL",
                combined,
            )

    def test_text_difference_fails(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = root / "expected"
            observed = root / "observed"

            self.write_bundle(
                expected,
                numeric="0.30",
                text="expected",
            )

            self.write_bundle(
                observed,
                numeric="0.30",
                text="changed",
            )

            result = self.run_comparison(
                expected,
                observed,
            )

            self.assertNotEqual(
                result.returncode,
                0,
            )

            combined = (
                result.stdout
                + result.stderr
            )

            self.assertIn(
                "column=label",
                combined,
            )


if __name__ == "__main__":
    unittest.main()
