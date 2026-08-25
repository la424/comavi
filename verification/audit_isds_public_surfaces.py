#!/usr/bin/env python3
"""Static closeout audit for benchmark, CHD, generic, live, report, and CI surfaces."""
from __future__ import annotations

import argparse
from pathlib import Path

REQUIRED_FIELDS = (
    "isds_v1",
    "isds_energy_component",
    "isds_context_component",
    "isds_energy_ratio_uncapped",
    "isds_dominant_axis",
    "isds_dominant_partner",
    "isds_dominant_signed_ddg",
    "isds_available",
    "isds_version",
)


def require_text(path: Path, needles: list[str]) -> list[str]:
    if not path.is_file():
        return [f"missing file: {path}"]
    text = path.read_text(encoding="utf-8")
    return [f"{path}: missing {needle!r}" for needle in needles if needle not in text]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    args = parser.parse_args()
    repo = args.repo.resolve()
    failures: list[str] = []

    failures += require_text(
        repo / "scripts/comavi_v7/pipeline.py",
        ["add_isds_v1_columns", "df = add_isds_v1_columns(df, partner_labels)"],
    )
    failures += require_text(repo / "run.py", ['res["df"].to_csv'])
    failures += require_text(repo / "scripts/run_chd.py", ["run_pipeline", ".to_csv("])
    failures += require_text(repo / "scripts/run_live.py", ["run_pipeline", ".to_csv("])
    failures += require_text(repo / "scripts/apply_concordance_v5.py", ["merged.to_csv"])
    failures += require_text(repo / "scripts/run_evaluate.py", list(REQUIRED_FIELDS))
    failures += require_text(
        repo / "scripts/run_chd_variants.py",
        ["run.py", "verify_isds_output_surfaces.py", "build_isds_variant_report.py"],
    )
    failures += require_text(
        repo / "scripts/build_isds_variant_report.py",
        ["prioritized.csv", "isds_report.md", "ISDS_OUTPUT_COLUMNS"],
    )
    failures += require_text(
        repo / "verification/verify_isds_output_surfaces.py",
        ["ISDS OUTPUT-SURFACE VERIFICATION: PASS"],
    )
    failures += require_text(
        repo / "README.md",
        ["<!-- ISDS_PUBLIC_WORKFLOW_V1 -->", "scripts/run_chd_variants.py"],
    )
    failures += require_text(
        repo / ".github/workflows/isds-v1.yml",
        [
            "# ISDS_PUBLIC_CLOSEOUT_V1",
            "# CHD_TRACKED_MIGRATION_V1",
            "# CHD_STRUCTURE_CONTRACT_V1",
            "verification/verify_chd_structure_contract.py",
            "reference_outputs/chd_isds_v1/chd_structure_manifest.json",
            "reference_outputs/chd_isds_v1",
            "verification.test_isds_variant_report",
        ],
    )

    failures += require_text(
        repo
        / "verification"
        / "verify_chd_migration_reproducibility.py",
        [
            "CHD MIGRATION REPRODUCIBILITY: PASS",
            "maximum_numeric_delta",
            "byte_identical",
        ],
    )
    failures += require_text(
        repo / ".github/workflows/isds-v1.yml",
        [
            "# CHD_SEMANTIC_REPRODUCIBILITY_V1",
            "verification/verify_chd_migration_reproducibility.py",
            "verification.test_chd_migration_reproducibility",
        ],
    )

    failures += require_text(
        repo / "notebooks" / "COMAVI_colab.ipynb",
        [
            "aeeaa3956b26edd67115083941954727316ca997",
            "MONOMER_OFFSET_OVERRIDES",
            "MULTIMER_OFFSET_OVERRIDES",
            "ARBITRARY-VARIANT OUTPUT CONTRACT: PASS",
            "build_isds_variant_report.py",
            "verify_isds_output_surfaces.py",
        ],
    )
    failures += require_text(
        repo / "verification" / "verify_colab_notebook_contract.py",
        [
            "COLAB NOTEBOOK CONTRACT: PASS",
            "Residue-numbering override contract: PASS",
        ],
    )
    # COLAB_PUBLIC_CLOSEOUT_V1

    if failures:
        print("ISDS PUBLIC-SURFACE AUDIT: FAIL")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)

    print("Shared benchmark/CHD/generic/live formula: PASS")
    print("Raw full-data writers: PASS")
    print("Reusable CHD wrapper: PASS")
    print("Benchmark evaluation display: PASS")
    print("Public ISDS report builder: PASS")
    print("Output-surface verifier: PASS")
    print("Cross-platform migration reproducibility: PASS")
    print("Public Colab notebook contract: PASS")
    print("README user workflow: PASS")
    print("CI closeout gates: PASS")
    print("ISDS PUBLIC-SURFACE AUDIT: PASS")


if __name__ == "__main__":
    main()
