#!/usr/bin/env python3
"""Verify cross-platform reproducibility of committed CHD ISDS-v1 reports.

CSV files are compared semantically because equivalent floating-point values
can be serialized with different final decimal digits across platforms.

The verifier still requires identical schema, column order, row order,
missingness, nonnumeric text, parsed JSON content, and Markdown content.
Numeric cells must agree within a fixed absolute tolerance.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


CSV_FILES = (
    "chd_concordance_with_isds_v1.csv",
    "chd_concordance_prioritized.csv",
)

JSON_FILE = "chd_concordance_isds_summary.json"
MARKDOWN_FILE = "chd_concordance_isds_report.md"

DEFAULT_ATOL = 1e-12


@dataclass(frozen=True)
class CSVComparison:
    rows: int
    columns: int
    maximum_numeric_delta: float
    byte_identical: bool


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def normalized_text(
    series: pd.Series,
) -> pd.Series:
    return (
        series.astype("string")
        .fillna("")
    )


def numeric_view(
    series: pd.Series,
) -> tuple[pd.Series, pd.Series, bool]:
    text = normalized_text(series).str.strip()
    present = text.ne("")

    numeric = pd.to_numeric(
        text.where(present),
        errors="coerce",
    )

    numeric_like = (
        bool(present.any())
        and bool(
            numeric.loc[present]
            .notna()
            .all()
        )
    )

    return numeric, present, numeric_like


def identity_label(
    frame: pd.DataFrame,
    row_index: int,
) -> str:
    fields: list[str] = []

    for column in (
        "system",
        "gene",
        "variant",
        "ref_aa",
        "position",
        "alt_aa",
    ):
        if column not in frame.columns:
            continue

        value = str(
            frame.iloc[row_index][column]
        ).strip()

        if not value:
            continue

        fields.append(
            f"{column}={value}"
        )

    if fields:
        return ", ".join(fields)

    # Add 2 because CSV line 1 is the header.
    return f"csv_line={row_index + 2}"


def compare_csv(
    expected: Path,
    observed: Path,
    *,
    atol: float,
) -> CSVComparison:
    left = pd.read_csv(
        expected,
        low_memory=False,
        keep_default_na=False,
    )

    right = pd.read_csv(
        observed,
        low_memory=False,
        keep_default_na=False,
    )

    if list(left.columns) != list(right.columns):
        raise ValueError(
            f"{expected.name}: column names or order differ.\n"
            f"expected={list(left.columns)}\n"
            f"observed={list(right.columns)}"
        )

    if len(left) != len(right):
        raise ValueError(
            f"{expected.name}: row count differs: "
            f"expected={len(left)}, observed={len(right)}"
        )

    failures: list[str] = []
    maximum_numeric_delta = 0.0

    for column in left.columns:
        (
            left_numeric,
            left_present,
            left_is_numeric,
        ) = numeric_view(left[column])

        (
            right_numeric,
            right_present,
            right_is_numeric,
        ) = numeric_view(right[column])

        if left_is_numeric and right_is_numeric:
            if not left_present.equals(
                right_present
            ):
                indexes = np.flatnonzero(
                    (
                        left_present
                        != right_present
                    ).to_numpy()
                )

                for index in indexes[:10]:
                    failures.append(
                        f"{identity_label(left, int(index))}; "
                        f"column={column}; "
                        "numeric missingness differs"
                    )

                continue

            left_values = left_numeric.to_numpy(
                dtype=float,
                na_value=np.nan,
            )

            right_values = right_numeric.to_numpy(
                dtype=float,
                na_value=np.nan,
            )

            match = np.isclose(
                left_values,
                right_values,
                atol=atol,
                rtol=0.0,
                equal_nan=True,
            )

            finite = (
                np.isfinite(left_values)
                & np.isfinite(right_values)
            )

            if finite.any():
                maximum_numeric_delta = max(
                    maximum_numeric_delta,
                    float(
                        np.max(
                            np.abs(
                                left_values[finite]
                                - right_values[finite]
                            )
                        )
                    ),
                )

            for index in np.flatnonzero(
                ~match
            )[:10]:
                failures.append(
                    f"{identity_label(left, int(index))}; "
                    f"column={column}; "
                    f"expected={left.iloc[index][column]!r}; "
                    f"observed={right.iloc[index][column]!r}; "
                    "delta="
                    f"{abs(left_values[index] - right_values[index]):.3e}"
                )

            continue

        left_text = normalized_text(
            left[column]
        )

        right_text = normalized_text(
            right[column]
        )

        mismatch = left_text.ne(
            right_text
        )

        for index in np.flatnonzero(
            mismatch.to_numpy()
        )[:10]:
            failures.append(
                f"{identity_label(left, int(index))}; "
                f"column={column}; "
                f"expected={left_text.iloc[index]!r}; "
                f"observed={right_text.iloc[index]!r}"
            )

    if failures:
        details = "\n  ".join(
            failures[:25]
        )

        raise ValueError(
            f"{expected.name}: semantic differences detected "
            "(showing up to 25):\n"
            f"  {details}"
        )

    return CSVComparison(
        rows=len(left),
        columns=len(left.columns),
        maximum_numeric_delta=(
            maximum_numeric_delta
        ),
        byte_identical=(
            sha256(expected)
            == sha256(observed)
        ),
    )


def normalize_newlines(
    text: str,
) -> str:
    return (
        text.replace("\r\n", "\n")
        .replace("\r", "\n")
    )


def compare_json(
    expected: Path,
    observed: Path,
) -> None:
    left = json.loads(
        expected.read_text(
            encoding="utf-8"
        )
    )

    right = json.loads(
        observed.read_text(
            encoding="utf-8"
        )
    )

    if left != right:
        raise ValueError(
            f"{expected.name}: parsed JSON objects differ"
        )


def compare_markdown(
    expected: Path,
    observed: Path,
) -> None:
    left = normalize_newlines(
        expected.read_text(
            encoding="utf-8"
        )
    )

    right = normalize_newlines(
        observed.read_text(
            encoding="utf-8"
        )
    )

    if left == right:
        return

    difference = "".join(
        difflib.unified_diff(
            left.splitlines(
                keepends=True
            ),
            right.splitlines(
                keepends=True
            ),
            fromfile=str(expected),
            tofile=str(observed),
            n=3,
        )
    )

    raise ValueError(
        f"{expected.name}: Markdown differs:\n"
        f"{difference[:12000]}"
    )


def require_files(
    root: Path,
    names: Iterable[str],
) -> None:
    missing = [
        name
        for name in names
        if not (root / name).is_file()
    ]

    if missing:
        raise ValueError(
            f"Missing files under {root}: "
            f"{missing}"
        )


def compare_directories(
    expected_dir: Path,
    observed_dir: Path,
    *,
    atol: float,
) -> None:
    expected_dir = (
        expected_dir.resolve()
    )

    observed_dir = (
        observed_dir.resolve()
    )

    required = (
        *CSV_FILES,
        JSON_FILE,
        MARKDOWN_FILE,
    )

    require_files(
        expected_dir,
        required,
    )

    require_files(
        observed_dir,
        required,
    )

    for name in CSV_FILES:
        result = compare_csv(
            expected_dir / name,
            observed_dir / name,
            atol=atol,
        )

        print(
            f"PASS {name}: "
            f"rows={result.rows}; "
            f"columns={result.columns}; "
            "maximum_numeric_delta="
            f"{result.maximum_numeric_delta:.3e}; "
            f"byte_identical={result.byte_identical}"
        )

    compare_json(
        expected_dir / JSON_FILE,
        observed_dir / JSON_FILE,
    )

    print(
        f"PASS {JSON_FILE}: "
        "parsed objects identical"
    )

    compare_markdown(
        expected_dir / MARKDOWN_FILE,
        observed_dir / MARKDOWN_FILE,
    )

    print(
        f"PASS {MARKDOWN_FILE}: "
        "normalized text identical"
    )

    print(
        "CHD MIGRATION REPRODUCIBILITY: PASS"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "--expected-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--observed-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--atol",
        type=float,
        default=DEFAULT_ATOL,
    )

    args = parser.parse_args()

    try:
        compare_directories(
            args.expected_dir,
            args.observed_dir,
            atol=args.atol,
        )
    except ValueError as error:
        raise SystemExit(
            "CHD MIGRATION REPRODUCIBILITY: FAIL\n"
            f"{error}"
        ) from error


if __name__ == "__main__":
    main()
