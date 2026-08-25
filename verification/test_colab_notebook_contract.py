from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "COMAVI_colab.ipynb"
VERIFIER = ROOT / "verification" / "verify_colab_notebook_contract.py"


class TestColabNotebookContract(unittest.TestCase):
    def run_verifier(self, notebook: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(VERIFIER), str(notebook)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

    def mutated_notebook(self, old: str, new: str, destination: Path) -> None:
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        changed = 0
        for cell in notebook["cells"]:
            value = cell.get("source", [])
            text = value if isinstance(value, str) else "".join(value)
            if old in text:
                text = text.replace(old, new)
                changed += 1
            cell["source"] = text.splitlines(keepends=True)
        self.assertGreater(changed, 0, f"mutation anchor absent: {old}")
        destination.write_text(json.dumps(notebook), encoding="utf-8")

    def test_current_notebook_passes(self) -> None:
        result = self.run_verifier(NOTEBOOK)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("COLAB NOTEBOOK CONTRACT: PASS", result.stdout)
        self.assertIn("Residue-numbering override contract: PASS", result.stdout)

    def test_stale_release_ref_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "stale.ipynb"
            self.mutated_notebook(
                "aeeaa3956b26edd67115083941954727316ca997",
                "v1.2-methods",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)

    def test_missing_isds_field_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "missing_isds.ipynb"
            self.mutated_notebook(
                "isds_energy_ratio_uncapped",
                "removed_energy_ratio_field",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)

    def test_missing_numbering_override_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "missing_offset.ipynb"
            self.mutated_notebook(
                "MONOMER_OFFSET_OVERRIDES",
                "REMOVED_MONOMER_OFFSET_FIELD",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
