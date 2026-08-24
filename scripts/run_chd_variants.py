#!/usr/bin/env python3
"""Run COMAVI on a CHD-focused variant batch and finalize ISDS-v1 outputs.

This is a thin public-facing wrapper around the generic ``run.py`` entry point.
It does not implement or copy the COMAVI formula. After structural scoring, it
verifies the nine versioned ISDS-v1 fields and writes a transparent report.

Use ``scripts/run_chd.py`` to reproduce the frozen paper cohort. Use this
wrapper for new CHD variants supplied by a user. For non-CHD systems, call
``run.py`` directly with the appropriate YAML configuration.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO / "configs" / "chd_systems.yaml"


@dataclass(frozen=True)
class RunnerOptions:
    variants: Path
    structures: Path
    output_dir: Path
    config: Path = DEFAULT_CONFIG
    foldx: str = os.environ.get("FOLDX_BINARY", "foldx")
    n_runs: int = 5
    no_fanout: bool = False
    dry_run: bool = False
    skip_numbering_check: bool = False
    quiet: bool = False
    no_report: bool = False
    report_prefix: str = "chd"
    top_n: int = 25


def build_pipeline_command(options: RunnerOptions) -> list[str]:
    command = [
        sys.executable,
        str(REPO / "run.py"),
        "--config",
        str(options.config),
        "--variants",
        str(options.variants),
        "--structures",
        str(options.structures),
        "--out",
        str(options.output_dir),
        "--foldx",
        str(options.foldx),
        "--n-runs",
        str(options.n_runs),
    ]
    if options.no_fanout:
        command.append("--no-fanout")
    if options.dry_run:
        command.append("--dry-run")
    if options.skip_numbering_check:
        command.append("--skip-numbering-check")
    if options.quiet:
        command.append("--quiet")
    return command


def build_postprocess_commands(options: RunnerOptions) -> list[list[str]]:
    result_csv = options.output_dir / "structural_results.csv"
    commands = [
        [
            sys.executable,
            str(REPO / "verification" / "verify_isds_output_surfaces.py"),
            "--require-available",
            str(result_csv),
        ]
    ]
    if not options.no_report:
        commands.append(
            [
                sys.executable,
                str(REPO / "scripts" / "build_isds_variant_report.py"),
                str(result_csv),
                "--out-dir",
                str(options.output_dir / "isds_v1_report"),
                "--prefix",
                options.report_prefix,
                "--top-n",
                str(options.top_n),
            ]
        )
    return commands


def parse_args() -> RunnerOptions:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, dest="output_dir")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--foldx", default=os.environ.get("FOLDX_BINARY", "foldx"))
    parser.add_argument("--n-runs", type=int, default=5)
    parser.add_argument("--no-fanout", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-numbering-check", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--no-report", action="store_true")
    parser.add_argument("--report-prefix", default="chd")
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()
    if args.n_runs < 1:
        parser.error("--n-runs must be at least 1")
    if args.top_n < 1:
        parser.error("--top-n must be at least 1")
    return RunnerOptions(**vars(args))


def run_checked(command: list[str]) -> None:
    """Run one delegated command and propagate its exit status cleanly."""
    print("+", " ".join(command), flush=True)

    try:
        completed = subprocess.run(
            command,
            cwd=REPO,
            check=False,
        )
    except OSError as error:
        raise SystemExit(
            f"Could not start delegated command: {error}"
        ) from None

    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


def main() -> None:
    options = parse_args()
    run_checked(build_pipeline_command(options))
    if options.dry_run:
        print("Dry run complete; no structural or report CSV was generated.")
        return
    result_csv = options.output_dir / "structural_results.csv"
    if not result_csv.is_file():
        raise SystemExit(f"Expected COMAVI result CSV was not created: {result_csv}")
    for command in build_postprocess_commands(options):
        run_checked(command)
    print(f"Verified structural results: {result_csv}")
    if not options.no_report:
        print(f"ISDS-v1 report directory: {options.output_dir / 'isds_v1_report'}")


if __name__ == "__main__":
    main()
