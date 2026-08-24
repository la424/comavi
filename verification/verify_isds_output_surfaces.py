#!/usr/bin/env python3
"""Verify benchmark, CHD, or live/external COMAVI CSVs against ISDS-v1."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from comavi_v7.isds import (  # noqa: E402
    ISDS_OUTPUT_COLUMNS,
    add_isds_v1_columns,
    infer_partner_labels,
)

NUMERIC = (
    "isds_v1",
    "isds_energy_component",
    "isds_context_component",
    "isds_energy_ratio_uncapped",
    "isds_dominant_signed_ddg",
)
TEXT = ("isds_dominant_axis", "isds_dominant_partner", "isds_version")


def bool_series(series: pd.Series) -> np.ndarray:
    return series.fillna(False).map(
        lambda value: str(value).strip().lower() in {"true", "1", "1.0", "yes"}
    ).to_numpy(bool)


def text_series(series: pd.Series) -> np.ndarray:
    return series.astype("string").fillna("").str.strip().to_numpy()


def verify(path: Path, tolerance: float, require_available: bool) -> None:
    df = pd.read_csv(path, low_memory=False)
    missing = [c for c in ISDS_OUTPUT_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"FAIL {path}: missing fields {missing}")

    versions = sorted(set(df["isds_version"].dropna().astype(str)) - {"ISDS-v1"})
    if versions:
        raise SystemExit(f"FAIL {path}: unexpected versions {versions}")

    labels = infer_partner_labels(df.columns)
    base = df.drop(columns=list(ISDS_OUTPUT_COLUMNS))
    recomputed = add_isds_v1_columns(base, labels)

    maximum_delta = 0.0
    for column in NUMERIC:
        observed = pd.to_numeric(df[column], errors="coerce").to_numpy(float)
        expected = pd.to_numeric(recomputed[column], errors="coerce").to_numpy(float)
        if not np.allclose(observed, expected, atol=tolerance, rtol=0.0, equal_nan=True):
            finite = np.isfinite(observed) & np.isfinite(expected)
            delta = float(np.max(np.abs(observed[finite] - expected[finite]))) if finite.any() else float("nan")
            raise SystemExit(f"FAIL {path}: {column} differs; maximum delta={delta}")
        finite = np.isfinite(observed) & np.isfinite(expected)
        if finite.any():
            maximum_delta = max(maximum_delta, float(np.max(np.abs(observed[finite] - expected[finite]))))

    for column in TEXT:
        if not np.array_equal(text_series(df[column]), text_series(recomputed[column])):
            raise SystemExit(f"FAIL {path}: text field {column} differs")

    observed_available = bool_series(df["isds_available"])
    expected_available = bool_series(recomputed["isds_available"])
    if not np.array_equal(observed_available, expected_available):
        raise SystemExit(f"FAIL {path}: isds_available differs")

    if require_available and int(observed_available.sum()) == 0:
        raise SystemExit(f"FAIL {path}: no rows have available ISDS-v1")

    for column in ("isds_v1", "isds_energy_component", "isds_context_component"):
        unavailable_values = pd.to_numeric(df.loc[~observed_available, column], errors="coerce")
        if unavailable_values.notna().any():
            raise SystemExit(f"FAIL {path}: unavailable rows contain {column}")

    # Formula-level cohort invariance on one available row.
    available_indices = np.flatnonzero(observed_available)
    if len(available_indices):
        index = int(available_indices[0])
        single = add_isds_v1_columns(base.iloc[[index]].copy(), labels).iloc[0]["isds_v1"]
        batch = recomputed.iloc[index]["isds_v1"]
        if single != batch:
            raise SystemExit(f"FAIL {path}: cohort invariance failed at row {index}")

    print(
        f"PASS {path}: rows={len(df)}; available={int(observed_available.sum())}; "
        f"partners={len(labels)}; maximum_numeric_delta={maximum_delta:.3e}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, nargs="+")
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--require-available", action="store_true")
    args = parser.parse_args()
    for path in args.csv:
        verify(path, args.tolerance, args.require_available)
    print("ISDS OUTPUT-SURFACE VERIFICATION: PASS")


if __name__ == "__main__":
    main()
