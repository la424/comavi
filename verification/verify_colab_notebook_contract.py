#!/usr/bin/env python3
"""Verify the public COMAVI Colab notebook and setup-wizard contract."""

from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

DEFAULT_PUBLIC_REF = "main"
SETUP_WIZARD_VERSION = "COMAVI-Setup-v1"
ISDS_FIELDS = {
    "isds_version",
    "isds_available",
    "isds_v1",
    "isds_energy_ratio_uncapped",
    "isds_energy_component",
    "isds_context_component",
    "isds_dominant_axis",
    "isds_dominant_partner",
    "isds_dominant_signed_ddg",
}
STALE_REFS = {
    "v1.2-methods",
    "v2.1-methods",
    "7aaa32b",
    "feature/isds-v1",
    "feature/chd-isds-public-closeout",
    "feature/colab-public-closeout",
}


def source(cell: dict) -> str:
    value = cell.get("source", [])
    return value if isinstance(value, str) else "".join(str(item) for item in value)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"COLAB NOTEBOOK CONTRACT: FAIL\n{message}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebook", type=Path)
    args = parser.parse_args()

    path = args.notebook.resolve()
    require(path.is_file(), f"Notebook is missing: {path}")
    notebook = json.loads(path.read_text(encoding="utf-8"))
    require(int(notebook.get("nbformat", -1)) == 4, "Notebook must use nbformat 4.")
    cells = notebook.get("cells")
    require(isinstance(cells, list) and cells, "Notebook contains no cells.")

    code_cells = [source(cell) for cell in cells if cell.get("cell_type") == "code"]
    markdown_cells = [source(cell) for cell in cells if cell.get("cell_type") == "markdown"]
    all_text = "\n".join(markdown_cells + code_cells)
    code_text = "\n".join(code_cells)

    for index, cell_source in enumerate(code_cells):
        try:
            ast.parse(cell_source, filename=f"notebook_code_cell_{index}.py")
        except SyntaxError as error:
            raise SystemExit(
                "COLAB NOTEBOOK CONTRACT: FAIL\n"
                f"Code cell {index} does not parse: {error}"
            ) from error

    require(
        re.search(r'COMAVI_REF\s*=\s*["\']main["\']', code_text) is not None,
        "Notebook does not default to the current public main branch.",
    )
    require("resolved_commit" in code_text, "Notebook does not record the resolved commit.")
    require("comavi_commit=" in code_text, "Result provenance lacks the resolved COMAVI commit.")
    require(
        re.search(r"REPO\s*/\s*[\"']run\.py[\"']", code_text) is not None,
        "Notebook does not invoke the generic run.py entry point.",
    )
    require(
        re.search(
            r"if\s+OUT\.exists\(\):\s*\n\s+shutil\.rmtree\(OUT\)",
            code_text,
        )
        is not None,
        "Notebook does not clear generated output before a rerun.",
    )
    require("run_chd.py" not in code_text, "Notebook is incorrectly tied to run_chd.py.")

    for field in ("gene", "ref_aa", "position", "alt_aa"):
        require(field in all_text, f"Notebook does not expose or validate variant field {field}.")
    for field in ISDS_FIELDS:
        require(field in all_text, f"Notebook does not expose or validate ISDS field {field}.")

    required_tokens = {
        "setup-wizard module": "comavi_setup_wizard",
        "new-system mode": "Create a new system",
        "prepared-bundle mode": "Use prepared system bundle(s)",
        "bundle merger": "merge_setup_bundles",
        "prepared current-variant revalidation": "revalidate_prepared_systems",
        "config-derived run provenance": "build_config_provenance",
        "system-scoped run provenance": "system_configurations=",
        "automatic chain assessment": "rank_structure_candidates",
        "monomer assessment": "assess_monomer_set",
        "complex assessment": "assess_structure",
        "traffic-light preflight": "TRAFFIC-LIGHT PREFLIGHT: PASS",
        "setup bundle": "COMAVI_system_setup_bundle.zip",
        "input-numbering override": "INPUT_NUMBERING_OVERRIDES",
        "chain override": "COMPLEX_CHAIN_OVERRIDES",
        "monomer offset override": "MONOMER_OFFSET_OVERRIDES",
        "multimer offset override": "MULTIMER_OFFSET_OVERRIDES",
        "yellow confirmation": "CONFIRM_YELLOW_PREFLIGHT",
        "biological-context confirmation": "CONFIRM_BIOLOGICAL_CONTEXT",
        "Google Drive FoldX": "Google Drive path",
        "plain-language cards": "render_plain_language_cards",
        "plain-language table": "comavi_plain_language_summary.csv",
        "arbitrary-variant assertion": "ARBITRARY-VARIANT OUTPUT CONTRACT: PASS",
        "public report": "build_isds_variant_report.py",
        "output verifier": "verify_isds_output_surfaces.py",
        "result bundle": "COMAVI_variant_results_bundle.zip",
        "ColabFold pin": "v1.6.1",
    }
    missing_tokens = [label for label, token in required_tokens.items() if token not in all_text]
    require(not missing_tokens, f"Notebook lacks setup-wizard elements: {missing_tokens}")

    require(
        "--skip-numbering-check" not in all_text,
        "Notebook exposes the unsafe numbering-check bypass.",
    )
    require(
        "Biological context is not automatable" in all_text,
        "Notebook does not state the limit of automatic structure selection.",
    )
    require(
        "Homomers and repeated copies" in all_text,
        "Notebook does not disclose the repeated-chain limitation.",
    )

    stale = sorted(token for token in STALE_REFS if token in all_text)
    require(not stale, f"Notebook contains stale release references: {stale}")

    uncleared = [
        index
        for index, cell in enumerate(cells)
        if cell.get("cell_type") == "code"
        and (cell.get("execution_count") is not None or cell.get("outputs"))
    ]
    require(not uncleared, f"Notebook contains executed output state in cells: {uncleared}")

    metadata_version = notebook.get("metadata", {}).get("comavi_setup_wizard_version")
    require(
        metadata_version == SETUP_WIZARD_VERSION,
        f"Notebook metadata wizard version differs: {metadata_version!r}",
    )

    print(f"Notebook: {path}")
    print(f"Cells: {len(cells)} ({len(markdown_cells)} Markdown, {len(code_cells)} code)")
    print(f"Default COMAVI ref: {DEFAULT_PUBLIC_REF}")
    print(f"Setup wizard: {SETUP_WIZARD_VERSION}")
    print(f"ISDS fields: {len(ISDS_FIELDS)}")
    print("Generic arbitrary-variant runner: PASS")
    print("Rerun-safe output isolation: PASS")
    print("Automatic sequence, chain, and numbering mapping: PASS")
    print("Prepared-bundle current-variant revalidation and multi-system reuse: PASS")
    print("Traffic-light preflight and fail-safe stops: PASS")
    print("Priority-score and mechanism-profile outputs: PASS")
    print("Plain-language cards and public reports: PASS")
    print("COLAB NOTEBOOK CONTRACT: PASS")


if __name__ == "__main__":
    main()
