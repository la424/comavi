#!/usr/bin/env python3
"""Generate the ISDS tier-comparator variant-call table from the canonical.

Why this script exists
----------------------
`data/analysis_inputs/isds_v1/tier_comparator_variant_calls.csv` supplies the
POSITIVE CLASS (`structural_ground_truth`) for every discrimination statistic
in the prioritization layer: ROC AUC, average precision, top-k recovery, the
system-cluster bootstrap, and the marginal tier screens.

Until v7.7 that file was hand-maintained and frozen at an old commit. When the
v7.6/v7.7 ledger adjudication moved three variants into the positive class the
file did not follow, so the published prioritization statistics were computed
against superseded labels while appearing to be reproducible. Byte-identical
output across canonical revisions was not evidence of robustness -- it was
evidence the layer was disconnected from the canonical.

Every column here is now derived from `scored_61var_canonical.csv` using the
shipped definitions (`apply_concordance_v5.discover_partners`,
`verify_tier_construction.tier_score` / `.to_tier`), so the file cannot drift
from the canonical again.

Usage
-----
  PYTHONPATH=scripts python scripts/build_tier_comparator_calls.py \
      --canonical reference_outputs/scored_61var_canonical.csv \
      --out data/analysis_inputs/isds_v1/tier_comparator_variant_calls.csv

  # verify the committed file matches what the canonical implies (CI / review):
  PYTHONPATH=scripts python scripts/build_tier_comparator_calls.py --check
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd

import apply_concordance_v5 as ac
from verify_tier_construction import (
    STRUCTURAL_CLASSES,
    tier_score,
    to_tier,
)

SILENT_CLASS = "structurally_silent"
REPO = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
DEFAULT_OUT = REPO / "data" / "analysis_inputs" / "isds_v1" / "tier_comparator_variant_calls.csv"

COLUMNS = [
    "system",
    "variant",
    "expected_mech_class",
    "structural_ground_truth",
    "reported_tier",
    "reconstructed_full_score",
    "no_interface_score",
    "interface_partner_count",
    "full_tier_strong",
    "no_interface_tier_strong",
    "interface_only_positive",
]


def build(canonical: pathlib.Path) -> pd.DataFrame:
    df = pd.read_csv(canonical)
    partners = ac.discover_partners(df)

    pop = df[
        df["comavi_tier"].notna()
        & df["expected_mech_class"].isin(STRUCTURAL_CLASSES + [SILENT_CLASS])
    ].copy()

    full = pop.apply(lambda r: tier_score(r, partners, True), axis=1)
    noif = pop.apply(lambda r: tier_score(r, partners, False), axis=1)

    pop["reconstructed_full_score"] = [s for s, _ in full]
    pop["interface_partner_count"] = [n for _, n in full]
    pop["no_interface_score"] = [s for s, _ in noif]

    # Fidelity gate: the reimplemented score must reproduce the shipped tier.
    # If this fails the component formula has drifted from mechanism.py and no
    # downstream statistic can be trusted.
    recomputed = [to_tier(s) for s in pop["reconstructed_full_score"]]
    bad = [
        (s, v, a, b)
        for s, v, a, b in zip(pop["system"], pop["variant"], recomputed, pop["comavi_tier"])
        if a != b
    ]
    if bad:
        raise SystemExit(
            "tier reimplementation does not reproduce shipped comavi_tier for "
            f"{len(bad)}/{len(pop)} rows, e.g. {bad[:3]}"
        )

    pop["structural_ground_truth"] = pop["expected_mech_class"].isin(STRUCTURAL_CLASSES)
    pop["reported_tier"] = pop["comavi_tier"].astype(str)
    pop["full_tier_strong"] = pop["reported_tier"].isin(ac.FOOTPRINT_TIERS)
    pop["no_interface_tier_strong"] = [
        to_tier(s) in ac.FOOTPRINT_TIERS for s in pop["no_interface_score"]
    ]
    pop["interface_only_positive"] = pop["interface_partner_count"] > 0

    out = pop[COLUMNS].sort_values(["system", "variant"]).reset_index(drop=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--canonical", type=pathlib.Path, default=DEFAULT_CANON)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--check",
        action="store_true",
        help="compare the committed file against the canonical instead of writing",
    )
    args = ap.parse_args()

    built = build(args.canonical.expanduser().resolve())

    if args.check:
        have = pd.read_csv(args.out).sort_values(["system", "variant"]).reset_index(drop=True)
        key = ["system", "variant"]
        if set(map(tuple, have[key].values)) != set(map(tuple, built[key].values)):
            raise SystemExit(
                f"cohort mismatch: committed n={len(have)} vs canonical-derived n={len(built)}"
            )
        drift = []
        for col in COLUMNS:
            a, b = have[col].astype(str).values, built[col].astype(str).values
            n = int((a != b).sum())
            if n:
                drift.append(f"{col}: {n} rows differ")
        if drift:
            raise SystemExit(
                "committed tier-comparator table is stale relative to the canonical:\n  "
                + "\n  ".join(drift)
            )
        print(f"OK tier_comparator_variant_calls.csv matches canonical (n={len(built)}, "
              f"positives={int(built.structural_ground_truth.sum())})")
        return

    args.out.parent.mkdir(parents=True, exist_ok=True)
    built.to_csv(args.out, index=False)
    print(
        f"wrote {args.out.relative_to(REPO)}  n={len(built)}  "
        f"positives={int(built.structural_ground_truth.sum())}  "
        f"strong_tier={int(built.full_tier_strong.sum())}"
    )


if __name__ == "__main__":
    sys.exit(main())
