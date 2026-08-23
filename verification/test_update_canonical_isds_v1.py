from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.comavi_v7.isds import ISDS_OUTPUT_COLUMNS
from scripts.update_canonical_isds_v1 import semantically_equal


REPO = Path(__file__).resolve().parents[1]
CANONICAL = (
    REPO
    / "reference_outputs"
    / "scored_61var_canonical.csv"
)
UPDATER = (
    REPO
    / "scripts"
    / "update_canonical_isds_v1.py"
)

OLD_NOTE = (
    "LABEL CORRECTION: split from C697F. Ollila 2006: 7.1% MMR, "
    "Lützen proteasome rescue confirms fold destabilization."
)


class CanonicalUpdaterTests(unittest.TestCase):
    def run_updater(
        self,
        source: Path,
        directory: Path,
        prefix: str,
    ) -> tuple[Path, Path, Path]:
        output = directory / f"{prefix}.csv"
        audit = directory / f"{prefix}_audit.csv"
        summary = directory / f"{prefix}_summary.json"

        result = subprocess.run(
            [
                sys.executable,
                str(UPDATER),
                "--input",
                str(source),
                "--output",
                str(output),
                "--audit",
                str(audit),
                "--summary",
                str(summary),
            ],
            cwd=REPO,
            capture_output=True,
            text=True,
        )

        self.assertEqual(
            result.returncode,
            0,
            msg=(
                "Updater failed.\nSTDOUT:\n"
                + result.stdout
                + "\nSTDERR:\n"
                + result.stderr
            ),
        )

        return output, audit, summary

    def test_blank_and_missing_values_are_equivalent(self):
        self.assertTrue(semantically_equal(np.nan, ""))
        self.assertTrue(semantically_equal(pd.NA, "   "))
        self.assertFalse(semantically_equal(np.nan, "monomer"))

    def test_current_canonical_is_idempotent(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)

            output, audit, summary_path = self.run_updater(
                CANONICAL,
                directory,
                "idempotent",
            )

            self.assertEqual(
                CANONICAL.read_bytes(),
                output.read_bytes(),
            )

            audit_table = pd.read_csv(audit, low_memory=False)
            self.assertEqual(len(audit_table), 0)

            summary = json.loads(
                summary_path.read_text(encoding="utf-8")
            )

            self.assertEqual(
                summary["initial_curation_state"],
                "already_updated",
            )
            self.assertEqual(
                summary["pre_existing_cells_changed"],
                0,
            )
            self.assertEqual(
                summary["isds_available_rows"],
                49,
            )

    def test_v21_style_a636p_record_migrates_to_canonical(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)

            current = pd.read_csv(
                CANONICAL,
                low_memory=False,
            )

            source = current.drop(
                columns=[
                    column
                    for column in ISDS_OUTPUT_COLUMNS
                    if column in current.columns
                ]
            ).copy()

            mask = (
                source["system"].astype(str).eq("msh2_msh6")
                & source["variant"].astype(str).eq("A636P")
            )

            self.assertEqual(int(mask.sum()), 1)

            source.loc[mask, "expected_ddg_monomer"] = "unknown"
            source.loc[mask, "expected_ddg_fold_complex"] = "unknown"
            source.loc[mask, "evidence_axes"] = "binding"
            source.loc[mask, "notes"] = OLD_NOTE
            source.loc[
                mask,
                "evaluation_note",
            ] = "high_axis_agreement_but_inconsistent_synthesis"

            for tag in ("t10", "t15", "t20", "t25", "tSAP"):
                source.loc[
                    mask,
                    f"structural_agreement_d_{tag}",
                ] = 2
                source.loc[
                    mask,
                    f"structural_agreement_{tag}",
                ] = 1.0
                source.loc[
                    mask,
                    f"directional_agreement_d_{tag}",
                ] = 2
                source.loc[
                    mask,
                    f"directional_agreement_{tag}",
                ] = 1.0

            old_style = directory / "old_style.csv"
            source.to_csv(old_style, index=False)

            output, audit, summary_path = self.run_updater(
                old_style,
                directory,
                "migrated",
            )

            self.assertEqual(
                CANONICAL.read_bytes(),
                output.read_bytes(),
            )

            audit_table = pd.read_csv(audit, low_memory=False)
            self.assertEqual(len(audit_table), 25)

            self.assertTrue(
                audit_table["system"]
                .astype(str)
                .eq("msh2_msh6")
                .all()
            )
            self.assertTrue(
                audit_table["variant"]
                .astype(str)
                .eq("A636P")
                .all()
            )

            summary = json.loads(
                summary_path.read_text(encoding="utf-8")
            )

            self.assertEqual(
                summary["initial_curation_state"],
                "pre_update",
            )
            self.assertEqual(
                summary["pre_existing_cells_changed"],
                25,
            )
            self.assertEqual(
                summary["isds_available_rows"],
                49,
            )


if __name__ == "__main__":
    unittest.main()
