#!/usr/bin/env python3
"""Verify the public COMAVI Colab notebook's static runtime contract."""

from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

CURRENT_PUBLIC_REF = "aeeaa3956b26edd67115083941954727316ca997"
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
    parser.add_argument("--expected-ref", default=CURRENT_PUBLIC_REF)
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

    require(args.expected_ref in all_text, f"Notebook does not pin expected public ref {args.expected_ref}.")
    require("COMAVI_REF" in code_text, "Notebook lacks COMAVI_REF.")
    require("resolved_commit" in code_text, "Notebook does not record the resolved commit.")
    require(
        re.search(r"REPO\s*/\s*[\"']run\.py[\"']", code_text) is not None,
        "Notebook does not invoke the generic run.py entry point.",
    )
    require("run_chd.py" not in code_text, "Notebook is incorrectly tied to run_chd.py.")

    for field in ("gene", "ref_aa", "position", "alt_aa"):
        require(field in all_text, f"Notebook does not expose or validate variant field {field}.")
    for field in ISDS_FIELDS:
        require(field in all_text, f"Notebook does not expose or validate ISDS field {field}.")

    require("build_isds_variant_report.py" in code_text, "Notebook does not build the public ISDS report.")
    require("verify_isds_output_surfaces.py" in code_text, "Notebook does not verify the public output contract.")
    require("ARBITRARY-VARIANT OUTPUT CONTRACT: PASS" in code_text, "Notebook lacks its variant-output assertion.")
    require("COMAVI_variant_results_bundle.zip" in code_text, "Notebook lacks a downloadable verified result bundle.")
    require("EXPERIMENTAL_CHAIN_MAP" in code_text, "Notebook lacks explicit experimental chain mapping.")
    require("COLABFOLD_RELEASE" in code_text and "v1.6.1" in code_text, "Notebook does not pin ColabFold v1.6.1.")

    for token in (
        "MONOMER_OFFSET_OVERRIDES",
        "MULTIMER_OFFSET_OVERRIDES",
        "NUMBERING OFFSET CONTRACT: PASS",
        "structure_position = submitted_position - offset",
        "monomer_offsets=",
        "multimer_offsets=",
        "gene_chain_map=",
    ):
        require(token in all_text, f"Notebook lacks numbering/provenance contract token: {token}")

    stale = sorted(token for token in STALE_REFS if token in all_text)
    require(not stale, f"Notebook contains stale release references: {stale}")

    uncleared = [
        index
        for index, cell in enumerate(cells)
        if cell.get("cell_type") == "code"
        and (cell.get("execution_count") is not None or cell.get("outputs"))
    ]
    require(not uncleared, f"Notebook contains executed output state in cells: {uncleared}")

    print(f"Notebook: {path}")
    print(f"Cells: {len(cells)} ({len(markdown_cells)} Markdown, {len(code_cells)} code)")
    print(f"Pinned COMAVI ref: {args.expected_ref}")
    print(f"ISDS fields: {len(ISDS_FIELDS)}")
    print("Generic runner: PASS")
    print("Variant input contract: PASS")
    print("Residue-numbering override contract: PASS")
    print("Priority-score and mechanism-profile outputs: PASS")
    print("Public report and output verifier: PASS")
    print("Experimental/PDB and ColabFold routes: STATIC PASS")
    print("COLAB NOTEBOOK CONTRACT: PASS")


if __name__ == "__main__":
    main()
