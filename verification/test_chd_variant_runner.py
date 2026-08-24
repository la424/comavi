from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import run_chd_variants as runner  # noqa: E402
from run_chd_variants import RunnerOptions, build_pipeline_command, build_postprocess_commands  # noqa: E402


class CHDVariantRunnerTests(unittest.TestCase):
    def options(self, **changes) -> RunnerOptions:
        base = dict(
            variants=Path("inputs/new_variants.csv"),
            structures=Path("structures"),
            output_dir=Path("results/new_chd_batch"),
            config=Path("configs/chd_systems.yaml"),
            foldx="/opt/foldx",
            n_runs=5,
        )
        base.update(changes)
        return RunnerOptions(**base)

    def test_wrapper_delegates_to_generic_runner(self):
        command = build_pipeline_command(self.options())
        self.assertEqual(Path(command[1]).name, "run.py")
        self.assertIn("configs/chd_systems.yaml", command)
        self.assertIn("inputs/new_variants.csv", command)
        self.assertIn("results/new_chd_batch", command)

    def test_postprocessing_verifies_and_reports_one_csv(self):
        commands = build_postprocess_commands(self.options())
        self.assertEqual(len(commands), 2)
        expected = "results/new_chd_batch/structural_results.csv"
        self.assertIn(expected, commands[0])
        self.assertIn(expected, commands[1])
        self.assertIn("--require-available", commands[0])

    def test_no_report_keeps_strict_verification(self):
        commands = build_postprocess_commands(self.options(no_report=True))
        self.assertEqual(len(commands), 1)
        self.assertIn("verify_isds_output_surfaces.py", commands[0][1])

    def test_child_exit_code_is_propagated_without_calledprocesserror(self):
        command = [
            sys.executable,
            "run.py",
        ]

        completed = subprocess.CompletedProcess(
            command,
            3,
        )

        with patch.object(
            runner.subprocess,
            "run",
            return_value=completed,
        ):
            with self.assertRaises(SystemExit) as raised:
                runner.run_checked(command)

        self.assertEqual(
            raised.exception.code,
            3,
        )


if __name__ == "__main__":
    unittest.main()
