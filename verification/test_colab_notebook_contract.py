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
        self.assertIn(
            "Prepared-bundle current-variant revalidation and multi-system reuse: PASS",
            result.stdout,
        )

    def test_stale_release_ref_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "stale.ipynb"
            self.mutated_notebook('COMAVI_REF = "main"', 'COMAVI_REF = "v1.2-methods"', candidate)
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

    def test_missing_setup_wizard_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "missing_wizard.ipynb"
            self.mutated_notebook(
                "comavi_setup_wizard",
                "removed_setup_wizard",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)

    def test_missing_prepared_bundle_mode_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "missing_bundle.ipynb"
            self.mutated_notebook(
                "Use prepared system bundle(s)",
                "Removed prepared mode",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)

    def test_missing_config_derived_provenance_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "missing_config_provenance.ipynb"
            self.mutated_notebook(
                "build_config_provenance",
                "removed_config_provenance",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)

    def test_missing_rerun_output_cleanup_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "missing_rerun_cleanup.ipynb"
            self.mutated_notebook(
                "shutil.rmtree(OUT)",
                "removed_rerun_cleanup(OUT)",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)

    def test_missing_biological_context_confirmation_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "missing_context_confirmation.ipynb"
            self.mutated_notebook(
                "CONFIRM_BIOLOGICAL_CONTEXT",
                "REMOVED_CONTEXT_CONFIRMATION",
                candidate,
            )
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)

    def test_numbering_bypass_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "unsafe.ipynb"
            notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
            notebook["cells"][0]["source"].append("\n--skip-numbering-check\n")
            candidate.write_text(json.dumps(notebook), encoding="utf-8")
            result = self.run_verifier(candidate)
            self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
