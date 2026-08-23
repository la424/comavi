#!/usr/bin/env python3
"""Verify that benchmark, CHD, or live/external COMAVI CSVs expose ISDS-v1."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REQUIRED = [
    "isds_v1",
    "isds_energy_component",
    "isds_context_component",
    "isds_energy_ratio_uncapped",
    "isds_dominant_axis",
    "isds_dominant_partner",
    "isds_dominant_signed_ddg",
    "isds_available",
    "isds_version",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, nargs="+")
    args = parser.parse_args()
    failures = []
    for path in args.csv:
        df = pd.read_csv(path, low_memory=False)
        missing = [c for c in REQUIRED if c not in df.columns]
        wrong_version = []
        if "isds_version" in df.columns:
            wrong_version = sorted(set(df["isds_version"].dropna().astype(str)) - {"ISDS-v1"})
        if missing or wrong_version:
            failures.append((path, missing, wrong_version))
            continue
        n_available = int(df["isds_available"].fillna(False).astype(bool).sum())
        print(f"PASS {path}: {n_available}/{len(df)} rows have ISDS-v1 available")
    if failures:
        for path, missing, versions in failures:
            print(f"FAIL {path}: missing={missing}; unexpected_versions={versions}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
