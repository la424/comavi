#!/usr/bin/env python3
"""Build transparent ISDS-v1 public-facing reports from a COMAVI result CSV.

The input may be a newly generated benchmark/CHD/live output or a historical
pre-ISDS result table. Existing COMAVI columns are preserved. The report adds
or independently verifies the nine versioned ISDS-v1 fields, writes a complete
augmented CSV, a compact prioritized CSV, a machine-readable summary, and a
human-readable Markdown report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_SCRIPTS = SCRIPT_DIR
if str(REPO_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(REPO_SCRIPTS))

from comavi_v7.isds import (  # noqa: E402
    ISDS_OUTPUT_COLUMNS,
    ISDS_VERSION,
    add_isds_v1_columns,
    infer_partner_labels,
)

NUMERIC_ISDS = (
    "isds_v1",
    "isds_energy_component",
    "isds_context_component",
    "isds_energy_ratio_uncapped",
    "isds_dominant_signed_ddg",
)
TEXT_ISDS = (
    "isds_dominant_axis",
    "isds_dominant_partner",
    "isds_version",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def normalize_bool(series: pd.Series) -> pd.Series:
    return series.fillna(False).map(
        lambda value: str(value).strip().lower() in {"true", "1", "1.0", "yes"}
    ).astype(bool)


def normalize_text(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("").str.strip()


def compare_stored_and_recomputed(
    stored: pd.DataFrame,
    recomputed: pd.DataFrame,
    *,
    tolerance: float = 1e-8,
) -> dict[str, object]:
    """Compare existing public ISDS fields with an independent recomputation."""
    missing_stored = [c for c in ISDS_OUTPUT_COLUMNS if c not in stored.columns]
    missing_recomputed = [c for c in ISDS_OUTPUT_COLUMNS if c not in recomputed.columns]
    if missing_stored or missing_recomputed:
        raise ValueError(
            f"ISDS schema mismatch: stored_missing={missing_stored}; "
            f"recomputed_missing={missing_recomputed}"
        )

    maximum_delta = 0.0
    failures: list[str] = []

    for column in NUMERIC_ISDS:
        left = pd.to_numeric(stored[column], errors="coerce")
        right = pd.to_numeric(recomputed[column], errors="coerce")
        match = (left.isna() & right.isna()) | ((left - right).abs() <= tolerance)
        if not bool(match.all()):
            failures.append(column)
        finite = left.notna() & right.notna()
        if bool(finite.any()):
            maximum_delta = max(maximum_delta, float((left[finite] - right[finite]).abs().max()))

    for column in TEXT_ISDS:
        if not normalize_text(stored[column]).equals(normalize_text(recomputed[column])):
            failures.append(column)

    if not normalize_bool(stored["isds_available"]).equals(
        normalize_bool(recomputed["isds_available"])
    ):
        failures.append("isds_available")

    return {
        "all_fields_match": not failures,
        "failing_fields": sorted(set(failures)),
        "maximum_numeric_delta": maximum_delta,
        "tolerance": tolerance,
    }


def direct_energy_columns(columns: Iterable[str]) -> list[str]:
    excluded = (
        "_sd",
        "_runs",
        "_confident",
        "_ci95_",
        "_distinguishable_",
        "_indistinguishable",
        "_vote_",
    )
    selected = []
    for column in columns:
        if column == "ddg_monomer" or column.startswith(("ddg_fold_", "ddg_binding_")):
            if any(token in column for token in excluded):
                continue
            selected.append(column)
    preferred = ["ddg_monomer"]
    return [c for c in preferred if c in selected] + sorted(c for c in selected if c not in preferred)


def prioritized_columns(df: pd.DataFrame) -> list[str]:
    identity = [
        "system",
        "gene",
        "variant",
        "ref_aa",
        "position",
        "alt_aa",
        "role",
        "phenotype",
    ]
    core = list(ISDS_OUTPUT_COLUMNS) + [
        "comavi_tier",
        "comavi_score",
        "evaluability",
        "evaluability_reason",
        "comavi_mechanism_tSAP",
        "comavi_mechanism_partner_tSAP",
        "comavi_mechanism_t25",
        "comavi_mechanism_partner_t25",
        "comavi_mechanism_t10",
        "comavi_mechanism_partner_t10",
        "comavi_mechanism",
        "comavi_mechanism_partner",
        "foldx_errors",
    ]
    ordered = identity + core + direct_energy_columns(df.columns)
    seen: set[str] = set()
    selected: list[str] = []
    for column in ordered:
        if column in df.columns and column not in seen:
            selected.append(column)
            seen.add(column)
    return selected


def build_outputs(
    input_csv: Path,
    output_dir: Path,
    *,
    prefix: str | None = None,
    replace_existing: bool = False,
    top_n: int = 25,
) -> dict[str, Path]:
    input_csv = input_csv.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = prefix or input_csv.stem

    source = pd.read_csv(input_csv, low_memory=False)
    if source.columns.duplicated().any():
        duplicates = sorted(set(source.columns[source.columns.duplicated()]))
        raise ValueError(f"Input contains duplicate columns: {duplicates}")

    partner_labels = infer_partner_labels(source.columns)
    if not partner_labels and "ddg_monomer" not in source.columns:
        raise ValueError("No monomer or partner energetic columns were found")

    existing = [c for c in ISDS_OUTPUT_COLUMNS if c in source.columns]
    partial = existing and len(existing) != len(ISDS_OUTPUT_COLUMNS)
    if partial:
        missing = [c for c in ISDS_OUTPUT_COLUMNS if c not in source.columns]
        raise ValueError(f"Input contains an incomplete ISDS schema; missing={missing}")

    source_without_isds = source.drop(columns=list(ISDS_OUTPUT_COLUMNS), errors="ignore")
    recomputed = add_isds_v1_columns(source_without_isds, partner_labels)

    audit: dict[str, object]
    if existing:
        audit = compare_stored_and_recomputed(source, recomputed)
        if not audit["all_fields_match"] and not replace_existing:
            raise ValueError(
                "Stored ISDS fields differ from independent recomputation: "
                f"{audit['failing_fields']}"
            )
        final = recomputed if replace_existing else source.copy()
    else:
        audit = {
            "all_fields_match": True,
            "failing_fields": [],
            "maximum_numeric_delta": 0.0,
            "tolerance": 1e-8,
            "note": "Input had no ISDS fields; fields were added from the shared formula.",
        }
        final = recomputed

    final_columns = list(source_without_isds.columns) + list(ISDS_OUTPUT_COLUMNS)
    final = final[final_columns]

    available = normalize_bool(final["isds_available"])
    unavailable = ~available
    for column in ("isds_v1", "isds_energy_component", "isds_context_component"):
        if final.loc[unavailable, column].notna().any():
            raise ValueError(f"Unavailable rows contain non-missing {column}")

    full_path = output_dir / f"{prefix}_with_isds_v1.csv"
    priority_path = output_dir / f"{prefix}_prioritized.csv"
    summary_path = output_dir / f"{prefix}_isds_summary.json"
    markdown_path = output_dir / f"{prefix}_isds_report.md"

    final.to_csv(full_path, index=False)

    sort_frame = final.copy()
    sort_frame["_isds_available_sort"] = normalize_bool(sort_frame["isds_available"]).astype(int)
    sort_frame["_isds_score_sort"] = pd.to_numeric(sort_frame["isds_v1"], errors="coerce").fillna(-math.inf)
    sort_frame["_isds_ratio_sort"] = pd.to_numeric(
        sort_frame["isds_energy_ratio_uncapped"], errors="coerce"
    ).fillna(-math.inf)
    priority = sort_frame.sort_values(
        ["_isds_available_sort", "_isds_score_sort", "_isds_ratio_sort"],
        ascending=[False, False, False],
        kind="mergesort",
    )
    priority = priority[prioritized_columns(final)]
    priority.to_csv(priority_path, index=False)

    system_count = int(final["system"].nunique()) if "system" in final.columns else None
    summary = {
        "input_file": input_csv.name,
        "input_sha256": sha256(input_csv),
        "isds_version": ISDS_VERSION,
        "rows": int(len(final)),
        "systems": system_count,
        "partner_labels": partner_labels,
        "isds_available_rows": int(available.sum()),
        "isds_unavailable_rows": int(unavailable.sum()),
        "existing_isds_schema": bool(existing),
        "stored_vs_recomputed_audit": audit,
        "outputs": {
            "full_csv": full_path.name,
            "prioritized_csv": priority_path.name,
            "markdown_report": markdown_path.name,
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    top = priority.head(top_n)
    report_lines = [
        f"# {prefix}: ISDS-v1 structural-disruption prioritization report",
        "",
        f"- Input rows: **{len(final)}**",
        f"- ISDS-v1 available: **{int(available.sum())}/{len(final)}**",
        f"- Systems: **{system_count if system_count is not None else 'not available'}**",
        f"- Formula version: **{ISDS_VERSION}**",
        "",
        "ISDS-v1 is a unitless structural-disruption prioritization index. It is not a "
        "pathogenicity probability or a validated binary decision rule. Use the signed mechanism "
        "profile and model-scope fields to choose follow-up experiments.",
        "",
        f"## Top {min(top_n, len(top))} rows",
        "",
    ]
    display = [
        c
        for c in (
            "system",
            "gene",
            "variant",
            "isds_v1",
            "isds_energy_component",
            "isds_context_component",
            "isds_dominant_axis",
            "isds_dominant_partner",
            "isds_dominant_signed_ddg",
            "comavi_tier",
            "comavi_mechanism_tSAP",
            "comavi_mechanism_t25",
        )
        if c in top.columns
    ]
    if display:
        visible = top[display].copy()
        report_lines.append("| " + " | ".join(display) + " |")
        report_lines.append("| " + " | ".join("---" for _ in display) + " |")
        for _, row in visible.iterrows():
            values = []
            for column in display:
                value = row[column]
                rendered = "" if pd.isna(value) else str(value)
                rendered = rendered.replace("|", "\\|").replace("\n", " ")
                values.append(rendered)
            report_lines.append("| " + " | ".join(values) + " |")
    else:
        report_lines.append("No displayable prioritization columns were found.")
    report_lines.extend(
        [
            "",
            "## Output files",
            "",
            f"- Complete augmented table: `{full_path.name}`",
            f"- Compact prioritized table: `{priority_path.name}`",
            f"- Machine-readable summary: `{summary_path.name}`",
        ]
    )
    markdown_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    return {
        "full_csv": full_path,
        "prioritized_csv": priority_path,
        "summary_json": summary_path,
        "markdown_report": markdown_path,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix")
    parser.add_argument(
        "--replace-existing",
        action="store_true",
        help="Replace stale existing ISDS fields instead of failing on a mismatch.",
    )
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    try:
        outputs = build_outputs(
            args.input_csv,
            args.out_dir,
            prefix=args.prefix,
            replace_existing=args.replace_existing,
            top_n=args.top_n,
        )
    except (
        ValueError,
        FileNotFoundError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ) as error:
        parser.error(str(error))

    for label, path in outputs.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
