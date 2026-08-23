#!/usr/bin/env python3
"""Build or verify the authoritative current COMAVI results summary."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import apply_concordance_v5 as concordance  # noqa: E402


CANONICAL = (
    REPO
    / "reference_outputs"
    / "scored_61var_canonical.csv"
)

LEDGER = (
    REPO
    / "reference_outputs"
    / "COMAVI_numbers_ledger.json"
)

STRESS_STATISTICS = (
    REPO
    / "reference_outputs"
    / "stress_tests"
    / "COMAVI_stress_verified_statistics.json"
)

STRESS_SUMMARY = (
    REPO
    / "reference_outputs"
    / "stress_tests"
    / "comavi_stress_tests.csv"
)

ISDS_SUMMARY = (
    REPO
    / "reference_outputs"
    / "isds_v1"
    / "ISDS_v1_summary.json"
)

DEFAULT_OUTPUT = (
    REPO
    / "docs"
    / "COMAVI_current_results.md"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Destination Markdown file.",
    )

    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "Exit nonzero when the existing output differs "
            "from regenerated content."
        ),
    )

    return parser.parse_args()


def stress_row(
    table: pd.DataFrame,
    name: str,
) -> pd.Series:
    rows = table.loc[
        table["test"].astype(str).eq(name)
    ]

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one stress-test row for {name!r}; "
            f"found {len(rows)}."
        )

    return rows.iloc[0]


def build_text() -> str:
    canonical = pd.read_csv(
        CANONICAL,
        low_memory=False,
    )

    ledger = json.loads(
        LEDGER.read_text(encoding="utf-8")
    )

    stress = json.loads(
        STRESS_STATISTICS.read_text(
            encoding="utf-8"
        )
    )

    stress_table = pd.read_csv(
        STRESS_SUMMARY,
        low_memory=False,
    )

    isds = json.loads(
        ISDS_SUMMARY.read_text(encoding="utf-8")
    )

    if len(canonical) != 61:
        raise RuntimeError(
            f"Canonical table contains {len(canonical)} rows."
        )

    partners = [
        partner
        for partner in concordance.discover_partners(
            canonical
        )
        if partner
        and "_ci95_" not in partner
        and "_distinguishable_" not in partner
    ]

    primary_mask = canonical[
        "mech_consistency_t25"
    ].isin(
        {
            "consistent",
            "partial",
            "inconsistent",
        }
    )

    axis_names = (
        "tier",
        "monomer",
        "fold",
        "binding",
    )

    def tally(mask: pd.Series) -> dict[str, list[int]]:
        totals = {
            axis: [0, 0]
            for axis in axis_names
        }

        for _, row in canonical.loc[mask].iterrows():
            result = (
                concordance.structural_agreement_by_axis(
                    row,
                    partners,
                    2.5,
                    2.5,
                    2.5,
                )
            )

            for axis, (
                numerator,
                denominator,
            ) in result.items():
                totals[axis][0] += int(numerator)
                totals[axis][1] += int(denominator)

        return totals

    all_axes = tally(
        pd.Series(
            True,
            index=canonical.index,
        )
    )

    primary_axes = tally(primary_mask)

    all_total = [
        sum(value[0] for value in all_axes.values()),
        sum(value[1] for value in all_axes.values()),
    ]

    primary_total = [
        sum(value[0] for value in primary_axes.values()),
        sum(value[1] for value in primary_axes.values()),
    ]

    energetic_primary = [
        sum(
            primary_axes[axis][0]
            for axis in (
                "monomer",
                "fold",
                "binding",
            )
        ),
        sum(
            primary_axes[axis][1]
            for axis in (
                "monomer",
                "fold",
                "binding",
            )
        ),
    ]

    if all_total != [99, 133]:
        raise RuntimeError(
            f"Unexpected all-row total: {all_total}"
        )

    if primary_total != [99, 132]:
        raise RuntimeError(
            f"Unexpected primary total: {primary_total}"
        )

    if energetic_primary != [65, 85]:
        raise RuntimeError(
            "Unexpected primary energetic-axis total: "
            f"{energetic_primary}"
        )

    grade_map = {
        "consistent": 1.0,
        "partial": 0.5,
        "inconsistent": 0.0,
    }

    mechanism_values = canonical[
        "mech_consistency_t25"
    ].map(grade_map).dropna()

    if (
        float(mechanism_values.sum()),
        len(mechanism_values),
    ) != (41.0, 57):
        raise RuntimeError(
            "Unexpected mechanism-consistency result."
        )

    observed_isds = isds["observed"]
    isds_population = isds["population"]

    positive_count = int(
        isds_population["positive"]
    )

    negative_count = int(
        isds_population["negative"]
    )

    top_k = {
        int(row["k"]): row
        for row in isds["top_k"]
    }

    for k in (10, 20, 25):
        if k not in top_k:
            raise RuntimeError(
                f"ISDS summary lacks top-{k}."
            )

    permutation_sa = stress_row(
        stress_table,
        "permutation null (SA)",
    )

    loso_sa = stress_row(
        stress_table,
        "leave-one-system-out (SA range)",
    )

    noise_sa = stress_row(
        stress_table,
        "replicate noise (SA)",
    )

    a636p_rows = canonical.loc[
        canonical["system"].astype(str).eq(
            "msh2_msh6"
        )
        & canonical["variant"].astype(str).eq(
            "A636P"
        )
    ]

    if len(a636p_rows) != 1:
        raise RuntimeError(
            "Expected one msh2_msh6 A636P row."
        )

    a636p = a636p_rows.iloc[0]

    tier_gradient = ledger["tier_gradient"]
    rho, p_value = ledger["tier_spearman"]

    cluster_sa_reportable = stress[
        "system_cluster_ci95"
    ]["sa"]["reportable_3"]

    cluster_mc_reportable = stress[
        "system_cluster_ci95"
    ]["mc"]["reportable_3"]

    variant_sa_reportable = stress[
        "variant_bootstrap_ci95"
    ]["sa"]["reportable_3"]

    variant_mc_reportable = stress[
        "variant_bootstrap_ci95"
    ]["mc"]["reportable_3"]

    text = f"""# COMAVI current results and denominator conventions

This file is generated by `scripts/build_current_results_summary.py` from the
committed canonical table and versioned reference outputs. It is the
authoritative current numerical summary. Historical v7 and v7.3 derivation
records remain in the repository but do not override this file.

## Canonical resource

- **61 variants** across **14 protein systems**.
- **57 variants** are gradeable for whole-variant mechanism-pattern agreement.
- **47 interaction variants** form the structural-prioritization population:
  **{positive_count}** committed modeled structural mechanisms and
  **{negative_count}** variants curated to lack a committed lesion on the
  modeled COMAVI axes.
- ISDS-v1 is available for
  **{int(canonical["isds_available"].fillna(False).astype(bool).sum())}/61**
  canonical rows.

MSH2 A636P is curated as neutral on the three modeled energetic axes, with
`expected_mech_class = {a636p["expected_mech_class"]}`. Its whole-variant
mechanism grade remains `{a636p["mech_consistency_t25"]}`, and its ISDS-v1
value is {float(a636p["isds_v1"]):.8f}.

## Mechanism localization

At the 2.5 kcal/mol reproducibility reference:

- weighted whole-variant mechanism-pattern agreement is
  **41/57 = {mechanism_values.mean():.4f}**;
- direction-aware agreement across the three energetic axes is
  **{energetic_primary[0]}/{energetic_primary[1]} =
  {energetic_primary[0] / energetic_primary[1]:.4f}**;
- monomer-fold agreement is
  **{primary_axes["monomer"][0]}/{primary_axes["monomer"][1]}**;
- complex-context agreement is
  **{primary_axes["fold"][0]}/{primary_axes["fold"][1]}**;
- binding agreement is
  **{primary_axes["binding"][0]}/{primary_axes["binding"][1]}**.

The four-output continuity aggregate, which adds the structural-context tier,
is:

- **{primary_total[0]}/{primary_total[1]} =
  {primary_total[0] / primary_total[1]:.4f}** on the primary 57-variant
  population;
- **{all_total[0]}/{all_total[1]} =
  {all_total[0] / all_total[1]:.4f}** when all 61 resource rows are retained.

The all-row decomposition is tier
{all_axes["tier"][0]}/{all_axes["tier"][1]}, monomer fold
{all_axes["monomer"][0]}/{all_axes["monomer"][1]}, complex context
{all_axes["fold"][0]}/{all_axes["fold"][1]}, and binding
{all_axes["binding"][0]}/{all_axes["binding"][1]}. The sole extra all-row
output is the monomer-fold output of BRCA1 R1699Q, which is retained in the
resource but excluded from whole-variant grading by its curated role.

## Mechanism uncertainty and robustness

The committed stress-test analysis uses 2,000 variant-bootstrap draws and
2,000 system-cluster draws. The system bootstrap resamples all 14 protein
systems with seed 8.

- mechanism-consistency variant-bootstrap 95% CI:
  **{variant_mc_reportable}**;
- mechanism-consistency system-cluster 95% CI:
  **{cluster_mc_reportable}**;
- structural-agreement variant-bootstrap 95% CI:
  **{variant_sa_reportable}**;
- structural-agreement system-cluster 95% CI:
  **{cluster_sa_reportable}**;
- structural-agreement leave-one-system-out range:
  **{loso_sa["null_mean"]}**;
- structural-agreement permutation null:
  observed **{float(permutation_sa["observed"]):.4f}** versus
  **{float(permutation_sa["null_mean"]):.4f} ±
  {float(permutation_sa["null_sd"]):.4f}**,
  **p = {float(permutation_sa["p_value"]):.4f}**;
- replicate-noise structural agreement:
  **{float(noise_sa["null_mean"]):.4f} ±
  {float(noise_sa["null_sd"]):.4f}**.

The complete summary and raw seeded draws are stored under
`reference_outputs/stress_tests/`.

## ISDS-v1 structural-disruption prioritization

On the {positive_count + negative_count}-variant prioritization population:

- ISDS-v1 ROC AUC:
  **{float(observed_isds["isds_roc_auc"]):.4f}**;
- ISDS-v1 average precision:
  **{float(observed_isds["isds_average_precision"]):.4f}**;
- Top 10:
  **{int(top_k[10]["structural_mechanisms_in_top_k"])}**
  structural-mechanism variants;
- Top 20:
  **{int(top_k[20]["structural_mechanisms_in_top_k"])}/{positive_count}**
  recovered;
- Top 25:
  **{int(top_k[25]["structural_mechanisms_in_top_k"])}/{positive_count}**
  recovered.

These are internal, system-aware benchmark results. ISDS-v1 is a unitless
prioritization index, not a pathogenicity probability or an externally
validated decision rule. No binary ISDS cutoff is established in this release.

## Structural-context tier

Among the 49 tier-carrying interaction variants, the observed pathogenic
fractions are:

- Tier 1:
  {tier_gradient["Tier 1"][0]}/{tier_gradient["Tier 1"][1]}
  ({tier_gradient["Tier 1"][2] * 100:.0f}%);
- Tier 2:
  {tier_gradient["Tier 2"][0]}/{tier_gradient["Tier 2"][1]}
  ({tier_gradient["Tier 2"][2] * 100:.0f}%);
- Tier 3:
  {tier_gradient["Tier 3"][0]}/{tier_gradient["Tier 3"][1]}
  ({tier_gradient["Tier 3"][2] * 100:.0f}%);
- Tier 4:
  {tier_gradient["Tier 4"][0]}/{tier_gradient["Tier 4"][1]}
  ({tier_gradient["Tier 4"][2] * 100:.0f}%).

The rank correlation is **ρ = {rho:.4f}, p = {p_value}**. The tier contains
no ΔΔG term and should be reported with its components because interface
membership materially contributes to its performance.

## Interpretation

COMAVI predicts modeled structural disruption rather than pathogenicity. A
positive energetic call supplies a testable monomer-fold, complex-context, or
binding hypothesis. A silent or scope-limited result means only that no lesion
was detected in the supplied static structural model at the stated operating
point.

## Reproducible sources

- Canonical table:
  `reference_outputs/scored_61var_canonical.csv`
- ISDS-v1 outputs:
  `reference_outputs/isds_v1/`
- Stress-test outputs:
  `reference_outputs/stress_tests/`
- Numerical ledger:
  `reference_outputs/COMAVI_numbers_ledger.json`
- Canonical verifier:
  `verification/verify_stage6.py`
- Stress-statistics verifier:
  `verification/verify_stress_statistics.py`
"""

    return text


def main() -> None:
    args = parse_args()
    output = args.output.expanduser().resolve()
    generated = build_text()

    if args.check:
        if not output.is_file():
            raise SystemExit(
                f"FAIL: current-results file is missing: {output}"
            )

        existing = output.read_text(
            encoding="utf-8"
        )

        if existing != generated:
            raise SystemExit(
                "FAIL: docs/COMAVI_current_results.md "
                "differs from regenerated content."
            )

        print(
            "CURRENT RESULTS DOCUMENT: "
            "UP TO DATE"
        )
        return

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output.write_text(
        generated,
        encoding="utf-8",
    )

    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
