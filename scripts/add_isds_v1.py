#!/usr/bin/env python3
"""Add ISDS-v1 columns to an existing COMAVI CSV without rerunning FoldX.

This supports historical benchmark, live, and CHD result files generated before
ISDS-v1 was integrated into the shared COMAVI engine.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from comavi_v7.isds import add_isds_v1_columns, infer_partner_labels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument(
        "--partner-label",
        action="append",
        dest="partner_labels",
        help="Partner label to use; repeat as needed. By default labels are inferred from DDG columns.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv, low_memory=False)
    labels = args.partner_labels or infer_partner_labels(df.columns)
    if not labels and "ddg_monomer" not in df.columns:
        parser.error("No monomer or partner energetic columns were found")
    out = add_isds_v1_columns(df, labels)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False)
    available = int(out["isds_available"].fillna(False).astype(bool).sum())
    print(f"Wrote {args.output_csv} with ISDS-v1 available for {available}/{len(out)} rows")


if __name__ == "__main__":
    main()
