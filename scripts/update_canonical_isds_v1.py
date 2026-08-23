#!/usr/bin/env python3
"""Apply the approved A636P curation update and add ISDS-v1.

This updater is intentionally narrow. It may change only:

1. four approved MSH2 A636P curation fields;
2. A636P structural- and directional-agreement denominators and rates that
   depend directly on those expectations;
3. A636P's evaluation note;
4. the nine new ISDS-v1 output columns.

Every other pre-existing canonical-table cell must remain semantically
unchanged. The script exits nonzero if that contract is violated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import apply_concordance_v5 as concordance  # noqa: E402
from comavi_v7.isds import (  # noqa: E402
    ISDS_OUTPUT_COLUMNS,
    ISDS_VERSION,
    add_isds_v1_columns,
)


A636P_SYSTEM = "msh2_msh6"
A636P_VARIANT = "A636P"

OLD_CURATION = {
    "expected_ddg_monomer": "unknown",
    "expected_ddg_fold_complex": "unknown",
    "expected_ddg_binding": "neutral",
    "evidence_axes": "binding",
    "notes": (
        "LABEL CORRECTION: split from C697F. Ollila 2006: 7.1% MMR, "
        "Lützen proteasome rescue confirms fold destabilization."
    ),
}

FINAL_NOTE = (
    "Variant-specific MSH2 studies report preserved expression/stability "
    "and preserved MSH2-MSH6 interaction, with impaired mismatch binding "
    "or ATP-dependent release. Monomer-fold, complex-context, and binding "
    "expectations are neutral; the mismatch-processing defect lies outside "
    "the three modeled COMAVI energetic axes."
)

FINAL_CURATION = {
    "expected_ddg_monomer": "neutral",
    "expected_ddg_fold_complex": "neutral",
    "expected_ddg_binding": "neutral",
    "evidence_axes": "functional",
    "notes": FINAL_NOTE,
}

THRESHOLD_TAGS = [
    tag for tag, _ in concordance.THRESHOLD_SPECS
]

APPROVED_EXISTING_CHANGES = {
    "expected_ddg_monomer",
    "expected_ddg_fold_complex",
    "evidence_axes",
    "notes",
    "evaluation_note",
}

for threshold_tag in THRESHOLD_TAGS:
    APPROVED_EXISTING_CHANGES.update(
        {
            f"structural_agreement_d_{threshold_tag}",
            f"structural_agreement_{threshold_tag}",
            f"directional_agreement_d_{threshold_tag}",
            f"directional_agreement_{threshold_tag}",
        }
    )

EXPECTED_FIRST_UPDATE_CHANGES = set(APPROVED_EXISTING_CHANGES)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply the approved A636P benchmark update and append "
            "cohort-independent ISDS-v1 fields."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Existing canonical 61-variant CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Destination for the updated canonical CSV.",
    )
    parser.add_argument(
        "--audit",
        type=Path,
        required=True,
        help="Destination for the cell-level change audit CSV.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        required=True,
        help="Destination for the machine-readable update summary JSON.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow replacement of existing output, audit, and summary files.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def is_missing(value: Any) -> bool:
    try:
        result = pd.isna(value)
    except (TypeError, ValueError):
        return False

    return bool(result) if isinstance(result, (bool, np.bool_)) else False


def is_blank_or_missing(value: Any) -> bool:
    """Treat CSV blank fields and pandas missing values as equivalent."""
    if is_missing(value):
        return True

    return isinstance(value, str) and not value.strip()


def semantically_equal(left: Any, right: Any) -> bool:
    left_blank = is_blank_or_missing(left)
    right_blank = is_blank_or_missing(right)

    if left_blank and right_blank:
        return True

    if left_blank or right_blank:
        return False

    if isinstance(left, (int, float, np.integer, np.floating)) and isinstance(
        right,
        (int, float, np.integer, np.floating),
    ):
        return math.isclose(
            float(left),
            float(right),
            rel_tol=0.0,
            abs_tol=1e-12,
        )

    return str(left) == str(right)


def normalize_axis_list(value: Any) -> tuple[str, ...]:
    if is_missing(value):
        return ()

    text = str(value).strip()

    if not text or text.lower() in {"nan", "none", "n/a"}:
        return ()

    return tuple(
        sorted(
            token.strip()
            for token in text.split(",")
            if token.strip()
        )
    )


def threshold_values(specification: Any) -> tuple[float, float, float]:
    if isinstance(specification, dict):
        return (
            float(specification["monomer"]),
            float(specification["fold"]),
            float(specification["binding"]),
        )

    value = float(specification)
    return value, value, value


def classify_curation_state(row: pd.Series) -> str:
    old_matches = all(
        semantically_equal(row.get(column), expected)
        for column, expected in OLD_CURATION.items()
    )
    final_matches = all(
        semantically_equal(row.get(column), expected)
        for column, expected in FINAL_CURATION.items()
    )

    if old_matches:
        return "pre_update"

    if final_matches:
        return "already_updated"

    observed = {
        column: row.get(column)
        for column in sorted(
            set(OLD_CURATION) | set(FINAL_CURATION)
        )
    }

    raise SystemExit(
        "STOP: A636P curation fields match neither the approved source state "
        f"nor the approved final state: {observed}"
    )


def compare_existing_cells(
    before: pd.DataFrame,
    after: pd.DataFrame,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []

    if len(before) != len(after):
        raise SystemExit(
            f"STOP: row count changed from {len(before)} to {len(after)}."
        )

    for column in before.columns:
        if column not in after.columns:
            raise SystemExit(
                f"STOP: pre-existing column was removed: {column}"
            )

        for index in before.index:
            old_value = before.at[index, column]
            new_value = after.at[index, column]

            if semantically_equal(old_value, new_value):
                continue

            records.append(
                {
                    "row_index": int(index),
                    "system": after.at[index, "system"],
                    "variant": after.at[index, "variant"],
                    "column": column,
                    "before": old_value,
                    "after": new_value,
                }
            )

    return pd.DataFrame(
        records,
        columns=[
            "row_index",
            "system",
            "variant",
            "column",
            "before",
            "after",
        ],
    )


def main() -> None:
    args = parse_args()

    input_path = args.input.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    audit_path = args.audit.expanduser().resolve()
    summary_path = args.summary.expanduser().resolve()

    if not input_path.is_file():
        raise SystemExit(f"Input canonical table not found: {input_path}")

    for destination in (output_path, audit_path, summary_path):
        if destination.exists() and not args.force:
            raise SystemExit(
                f"STOP: destination already exists: {destination}"
            )
        destination.parent.mkdir(parents=True, exist_ok=True)

    source = pd.read_csv(input_path, low_memory=False)
    original_columns = list(source.columns)

    if len(source) != 61:
        raise SystemExit(
            f"STOP: expected 61 canonical rows; found {len(source)}."
        )

    keys = ["system", "variant"]

    if source.duplicated(keys).any():
        duplicates = source.loc[
            source.duplicated(keys, keep=False),
            keys,
        ]
        raise SystemExit(
            "STOP: duplicate system/variant keys:\n"
            + duplicates.to_string(index=False)
        )

    mask = (
        source["system"].astype(str).eq(A636P_SYSTEM)
        & source["variant"].astype(str).eq(A636P_VARIANT)
    )

    if int(mask.sum()) != 1:
        raise SystemExit(
            "STOP: expected one msh2_msh6 A636P row; "
            f"found {int(mask.sum())}."
        )

    index = source.index[mask][0]
    candidate = source.copy()
    initial_state = classify_curation_state(candidate.loc[index])

    for column, value in FINAL_CURATION.items():
        if column not in candidate.columns:
            raise SystemExit(
                f"STOP: required curation column is absent: {column}"
            )
        candidate.at[index, column] = value

    derived_class = concordance.derive_expected_mech_class(
        candidate.loc[index]
    )

    if derived_class != "structurally_silent":
        raise SystemExit(
            "STOP: corrected A636P expectations did not derive the expected "
            f"'structurally_silent' class: {derived_class!r}"
        )

    stored_class = candidate.at[index, "expected_mech_class"]

    if not semantically_equal(stored_class, derived_class):
        raise SystemExit(
            "STOP: A636P expected_mech_class would change unexpectedly: "
            f"stored={stored_class!r}, derived={derived_class!r}"
        )

    partners = [
        partner
        for partner in concordance.discover_partners(candidate)
        if "_ci95_" not in partner
        and "_distinguishable_" not in partner
    ]

    axis_status = concordance.classify_axis_status(
        candidate.loc[index]
    )

    for tag, specification in concordance.THRESHOLD_SPECS:
        monomer_threshold, fold_threshold, binding_threshold = (
            threshold_values(specification)
        )

        row = candidate.loc[index]

        structural_n, structural_d = (
            concordance.compute_structural_agreement(
                row,
                partners,
                monomer_threshold,
                fold_threshold,
                binding_threshold,
            )
        )

        stored_structural_n = row.get(
            f"structural_agreement_n_{tag}"
        )

        if not semantically_equal(
            stored_structural_n,
            structural_n,
        ):
            raise SystemExit(
                "STOP: the A636P structural-agreement numerator changed "
                f"unexpectedly at {tag}: "
                f"stored={stored_structural_n!r}, "
                f"derived={structural_n!r}"
            )

        candidate.at[
            index,
            f"structural_agreement_d_{tag}",
        ] = structural_d

        candidate.at[
            index,
            f"structural_agreement_{tag}",
        ] = (
            structural_n / structural_d
            if structural_d
            else np.nan
        )

        directional_full, directional_half, directional_d = (
            concordance.compute_directional_agreement(
                row,
                partners,
                monomer_threshold,
                fold_threshold,
                binding_threshold,
            )
        )

        stored_directional_full = row.get(
            f"directional_agreement_full_{tag}"
        )
        stored_directional_half = row.get(
            f"directional_agreement_half_{tag}"
        )

        if not semantically_equal(
            stored_directional_full,
            directional_full,
        ):
            raise SystemExit(
                "STOP: the A636P directional-agreement full-credit count "
                f"changed unexpectedly at {tag}: "
                f"stored={stored_directional_full!r}, "
                f"derived={directional_full!r}"
            )

        if not semantically_equal(
            stored_directional_half,
            directional_half,
        ):
            raise SystemExit(
                "STOP: the A636P directional-agreement half-credit count "
                f"changed unexpectedly at {tag}: "
                f"stored={stored_directional_half!r}, "
                f"derived={directional_half!r}"
            )

        candidate.at[
            index,
            f"directional_agreement_d_{tag}",
        ] = directional_d

        candidate.at[
            index,
            f"directional_agreement_{tag}",
        ] = (
            (directional_full + 0.5 * directional_half)
            / directional_d
            if directional_d
            else np.nan
        )

        mechanism_column = f"mech_{tag}"
        grade_column = f"mech_consistency_{tag}"
        false_positive_column = f"mech_false_positive_axes_{tag}"
        missed_positive_column = f"mech_missed_positive_axes_{tag}"

        grade, false_positive, missed_positive = (
            concordance.grade_mechanism_consistency(
                candidate.loc[index],
                candidate.at[index, mechanism_column],
                derived_class,
                axis_status,
            )
        )

        if not semantically_equal(
            candidate.at[index, grade_column],
            grade,
        ):
            raise SystemExit(
                f"STOP: A636P {grade_column} would change unexpectedly: "
                f"stored={candidate.at[index, grade_column]!r}, "
                f"derived={grade!r}"
            )

        if normalize_axis_list(
            candidate.at[index, false_positive_column]
        ) != tuple(sorted(false_positive or [])):
            raise SystemExit(
                f"STOP: A636P {false_positive_column} would change "
                "unexpectedly."
            )

        if normalize_axis_list(
            candidate.at[index, missed_positive_column]
        ) != tuple(sorted(missed_positive or [])):
            raise SystemExit(
                f"STOP: A636P {missed_positive_column} would change "
                "unexpectedly."
            )

        nbhd_mechanism_column = f"nbhd_mech_{tag}"
        nbhd_grade_column = f"nbhd_mech_consistency_{tag}"
        nbhd_false_column = (
            f"nbhd_mech_false_positive_axes_{tag}"
        )
        nbhd_missed_column = (
            f"nbhd_mech_missed_positive_axes_{tag}"
        )

        if nbhd_mechanism_column in candidate.columns:
            nbhd_grade, nbhd_false, nbhd_missed = (
                concordance.grade_mechanism_consistency(
                    candidate.loc[index],
                    candidate.at[index, nbhd_mechanism_column],
                    derived_class,
                    axis_status,
                )
            )

            if not semantically_equal(
                candidate.at[index, nbhd_grade_column],
                nbhd_grade,
            ):
                raise SystemExit(
                    f"STOP: A636P {nbhd_grade_column} would change "
                    "unexpectedly."
                )

            if normalize_axis_list(
                candidate.at[index, nbhd_false_column]
            ) != tuple(sorted(nbhd_false or [])):
                raise SystemExit(
                    f"STOP: A636P {nbhd_false_column} would change "
                    "unexpectedly."
                )

            if normalize_axis_list(
                candidate.at[index, nbhd_missed_column]
            ) != tuple(sorted(nbhd_missed or [])):
                raise SystemExit(
                    f"STOP: A636P {nbhd_missed_column} would change "
                    "unexpectedly."
                )

    if not hasattr(concordance, "compute_evaluation_note"):
        raise SystemExit(
            "STOP: apply_concordance_v5 lacks compute_evaluation_note."
        )

    candidate.at[index, "evaluation_note"] = (
        concordance.compute_evaluation_note(
            candidate.loc[index],
            (
                candidate.at[
                    index,
                    "structural_agreement_n_t25",
                ],
                candidate.at[
                    index,
                    "structural_agreement_d_t25",
                ],
            ),
            candidate.at[index, "mech_consistency_t25"],
        )
    )

    # Remove any pre-existing ISDS columns and regenerate them from source
    # fields, so re-running the updater never trusts stale score values.
    candidate = candidate.drop(
        columns=[
            column
            for column in ISDS_OUTPUT_COLUMNS
            if column in candidate.columns
        ],
        errors="ignore",
    )

    candidate = add_isds_v1_columns(candidate, partners)

    existing_changes = compare_existing_cells(
        source,
        candidate,
    )

    if not existing_changes.empty:
        outside_a636p = existing_changes[
            ~(
                existing_changes["system"]
                .astype(str)
                .eq(A636P_SYSTEM)
                & existing_changes["variant"]
                .astype(str)
                .eq(A636P_VARIANT)
            )
        ]

        if not outside_a636p.empty:
            raise SystemExit(
                "STOP: pre-existing cells outside A636P changed:\n"
                + outside_a636p.to_string(index=False)
            )

        unapproved = existing_changes[
            ~existing_changes["column"].isin(
                APPROVED_EXISTING_CHANGES
            )
        ]

        if not unapproved.empty:
            raise SystemExit(
                "STOP: unapproved A636P fields changed:\n"
                + unapproved.to_string(index=False)
            )

    changed_columns = set(existing_changes["column"])

    if initial_state == "pre_update":
        missing_expected = (
            EXPECTED_FIRST_UPDATE_CHANGES - changed_columns
        )
        unexpected_expected = (
            changed_columns - EXPECTED_FIRST_UPDATE_CHANGES
        )

        if missing_expected or unexpected_expected:
            raise SystemExit(
                "STOP: first-update change set differs from the approved "
                f"contract; missing={sorted(missing_expected)}, "
                f"unexpected={sorted(unexpected_expected)}"
            )

    if initial_state == "already_updated" and changed_columns:
        raise SystemExit(
            "STOP: an already-updated canonical table still required "
            f"pre-existing changes: {sorted(changed_columns)}"
        )

    missing_isds = [
        column
        for column in ISDS_OUTPUT_COLUMNS
        if column not in candidate.columns
    ]

    if missing_isds:
        raise SystemExit(
            f"STOP: ISDS output columns are missing: {missing_isds}"
        )

    available = (
        candidate["isds_available"]
        .fillna(False)
        .astype(bool)
    )

    if int(available.sum()) != 49:
        raise SystemExit(
            "STOP: expected ISDS-v1 availability for 49/61 rows; "
            f"observed {int(available.sum())}/61."
        )

    final_columns = original_columns + [
        column
        for column in ISDS_OUTPUT_COLUMNS
        if column not in original_columns
    ]

    candidate = candidate[final_columns]

    existing_changes.to_csv(audit_path, index=False)
    candidate.to_csv(output_path, index=False)

    summary = {
        "status": "pass",
        "isds_version": ISDS_VERSION,
        "initial_curation_state": initial_state,
        "input": str(input_path),
        "output": str(output_path),
        "input_sha256": sha256(input_path),
        "output_sha256": sha256(output_path),
        "rows": len(candidate),
        "pre_existing_cells_changed": len(existing_changes),
        "changed_existing_columns": sorted(changed_columns),
        "changes_outside_a636p": 0,
        "isds_fields_added": list(ISDS_OUTPUT_COLUMNS),
        "isds_available_rows": int(available.sum()),
        "isds_unavailable_rows": int((~available).sum()),
        "a636p_expected_mech_class": derived_class,
        "a636p_mech_consistency_t25": candidate.at[
            index,
            "mech_consistency_t25",
        ],
        "a636p_structural_agreement_t25": {
            "numerator": int(
                candidate.at[
                    index,
                    "structural_agreement_n_t25",
                ]
            ),
            "denominator": int(
                candidate.at[
                    index,
                    "structural_agreement_d_t25",
                ]
            ),
        },
        "a636p_isds_v1": float(
            candidate.at[index, "isds_v1"]
        ),
    }

    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print(f"Input:   {input_path}")
    print(f"Output:  {output_path}")
    print(f"Audit:   {audit_path}")
    print(f"Summary: {summary_path}")
    print(f"Initial curation state: {initial_state}")
    print(
        "Pre-existing cells changed: "
        f"{len(existing_changes)}"
    )
    print(
        "Changed existing columns: "
        + ", ".join(sorted(changed_columns))
    )
    print(
        "ISDS-v1 available: "
        f"{int(available.sum())}/{len(candidate)}"
    )
    print(
        "A636P structural agreement at 2.5: "
        f"{int(candidate.at[index, 'structural_agreement_n_t25'])}/"
        f"{int(candidate.at[index, 'structural_agreement_d_t25'])}"
    )
    print(
        "A636P mechanism grade at 2.5: "
        f"{candidate.at[index, 'mech_consistency_t25']}"
    )
    print(
        "A636P ISDS-v1: "
        f"{float(candidate.at[index, 'isds_v1']):.8f}"
    )
    print("MINIMAL CANONICAL UPDATE: PASS")


if __name__ == "__main__":
    main()
