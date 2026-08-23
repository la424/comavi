#!/usr/bin/env python3
"""Verify the committed COMAVI stress-test statistics and provenance."""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]

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

STRESS_SCRIPT = (
    REPO
    / "verification"
    / "stress_tests.py"
)

ROOT = (
    REPO
    / "reference_outputs"
    / "stress_tests"
)

STATISTICS = (
    ROOT
    / "COMAVI_stress_verified_statistics.json"
)

SUMMARY = (
    ROOT
    / "comavi_stress_tests.csv"
)

DRAWS = (
    ROOT
    / "comavi_stress_draws.npz"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def require(
    condition: bool,
    message: str,
) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def close(
    left: float,
    right: float,
    tolerance: float = 1e-12,
) -> bool:
    return math.isclose(
        float(left),
        float(right),
        rel_tol=0.0,
        abs_tol=tolerance,
    )


for required_path in (
    CANONICAL,
    LEDGER,
    STRESS_SCRIPT,
    STATISTICS,
    SUMMARY,
    DRAWS,
):
    require(
        required_path.is_file(),
        f"required file is missing: {required_path}",
    )


statistics = json.loads(
    STATISTICS.read_text(encoding="utf-8")
)

ledger = json.loads(
    LEDGER.read_text(encoding="utf-8")
)

canonical = pd.read_csv(
    CANONICAL,
    low_memory=False,
)

summary = pd.read_csv(
    SUMMARY,
    low_memory=False,
)


# ---------------------------------------------------------------------
# Population and point estimates
# ---------------------------------------------------------------------

population = statistics["population"]
observed = statistics["observed"]
draw_metadata = statistics["draws"]

require(
    len(canonical) == 61,
    f"canonical table has {len(canonical)} rows, not 61",
)

require(
    int(population["canonical_rows"]) == 61,
    "statistics record canonical-row count differs from 61",
)

require(
    int(population["mechanism_graded_variants"]) == 57,
    "mechanism population differs from 57",
)

require(
    list(population["primary_structural_agreement"])
    == [99, 132],
    "primary structural-agreement count differs from 99/132",
)

require(
    list(population["all_row_structural_agreement"])
    == [99, 133],
    "all-row structural-agreement count differs from 99/133",
)

require(
    int(population["systems"]) == 14,
    "system count differs from 14",
)

require(
    close(
        observed["mechanism_consistency"],
        41 / 57,
    ),
    "mechanism-consistency point estimate differs from 41/57",
)

require(
    close(
        observed["structural_agreement"],
        99 / 132,
    ),
    "structural-agreement point estimate differs from 99/132",
)


# ---------------------------------------------------------------------
# File provenance
# ---------------------------------------------------------------------

expected_hashes = {
    "canonical_sha256": sha256(CANONICAL),
    "stress_script_sha256": sha256(STRESS_SCRIPT),
    "stress_summary_sha256": sha256(SUMMARY),
    "stress_draws_sha256": sha256(DRAWS),
}

for field, expected in expected_hashes.items():
    observed_hash = statistics["files"].get(field)

    require(
        observed_hash == expected,
        (
            f"{field} differs: "
            f"record={observed_hash!r}, current={expected!r}"
        ),
    )


# ---------------------------------------------------------------------
# Raw draws and interval reconstruction
# ---------------------------------------------------------------------

with np.load(DRAWS) as archive:
    required_keys = {
        "boot_mc",
        "boot_sa",
        "cluster_mc",
        "cluster_sa",
        "seed_cluster",
        "cluster_system_count",
        "cluster_systems",
        "mc_n",
        "sa_numerator",
        "sa_denominator",
        "mc_obs",
        "sa_obs",
    }

    missing = required_keys - set(archive.files)

    require(
        not missing,
        f"draw archive is missing fields: {sorted(missing)}",
    )

    arrays = {
        "variant_mc": archive["boot_mc"].astype(float),
        "variant_sa": archive["boot_sa"].astype(float),
        "cluster_mc": archive["cluster_mc"].astype(float),
        "cluster_sa": archive["cluster_sa"].astype(float),
    }

    seed_cluster = int(
        archive["seed_cluster"][0]
    )

    system_count = int(
        archive["cluster_system_count"][0]
    )

    system_labels = (
        archive["cluster_systems"]
        .astype(str)
        .tolist()
    )

    archive_mc_n = int(
        archive["mc_n"][0]
    )

    archive_sa_n = int(
        archive["sa_numerator"][0]
    )

    archive_sa_d = int(
        archive["sa_denominator"][0]
    )

    archive_mc = float(
        archive["mc_obs"][0]
    )

    archive_sa = float(
        archive["sa_obs"][0]
    )


require(
    seed_cluster == 8,
    f"cluster seed is {seed_cluster}, not 8",
)

require(
    system_count == 14,
    f"draw archive system count is {system_count}, not 14",
)

require(
    archive_mc_n == 57,
    f"draw archive mechanism n is {archive_mc_n}, not 57",
)

require(
    (archive_sa_n, archive_sa_d) == (99, 132),
    (
        "draw archive structural agreement is "
        f"{archive_sa_n}/{archive_sa_d}, not 99/132"
    ),
)

require(
    close(archive_mc, 41 / 57),
    "draw archive mechanism point estimate differs from 41/57",
)

require(
    close(archive_sa, 99 / 132),
    "draw archive structural point estimate differs from 99/132",
)

require(
    len(arrays["variant_sa"])
    == int(draw_metadata["variant_bootstrap"]),
    "variant-bootstrap draw count differs from the statistics record",
)

require(
    len(arrays["cluster_sa"])
    == int(draw_metadata["system_cluster_bootstrap"]),
    "cluster-bootstrap draw count differs from the statistics record",
)

require(
    int(draw_metadata["cluster_seed"]) == seed_cluster,
    "statistics-record cluster seed differs from the draw archive",
)

require(
    system_labels
    == list(population["system_labels"]),
    "system-label order differs between statistics record and draw archive",
)

require(
    set(system_labels)
    == set(canonical["system"].astype(str).unique()),
    "draw-archive systems differ from the canonical table",
)


def summary_row(test_name: str) -> pd.Series:
    rows = summary.loc[
        summary["test"].astype(str).eq(test_name)
    ]

    require(
        len(rows) == 1,
        (
            f"expected one summary row for {test_name!r}; "
            f"found {len(rows)}"
        ),
    )

    return rows.iloc[0]


interval_contracts: dict[
    str,
    tuple[str, str, str],
] = {
    "variant_mc": (
        "variant_bootstrap_ci95",
        "mc",
        "bootstrap 95% CI (MC)",
    ),
    "variant_sa": (
        "variant_bootstrap_ci95",
        "sa",
        "bootstrap 95% CI (SA)",
    ),
    "cluster_mc": (
        "system_cluster_ci95",
        "mc",
        "cluster bootstrap 95% CI (MC)",
    ),
    "cluster_sa": (
        "system_cluster_ci95",
        "sa",
        "cluster bootstrap 95% CI (SA)",
    ),
}


for name, (
    section,
    metric,
    test_name,
) in interval_contracts.items():
    calculated = np.nanpercentile(
        arrays[name],
        [2.5, 97.5],
    )

    recorded = statistics[section][metric]

    full_precision = np.asarray(
        recorded["full_precision"],
        dtype=float,
    )

    require(
        np.allclose(
            calculated,
            full_precision,
            atol=1e-15,
            rtol=0.0,
            equal_nan=True,
        ),
        (
            f"{name} full-precision interval differs: "
            f"calculated={calculated.tolist()}, "
            f"recorded={full_precision.tolist()}"
        ),
    )

    rounded_four = [
        round(float(value), 4)
        for value in calculated
    ]

    require(
        rounded_four
        == list(recorded["rounded_4"]),
        (
            f"{name} four-decimal interval differs: "
            f"{rounded_four} versus {recorded['rounded_4']}"
        ),
    )

    reportable_three = [
        round(float(value), 3)
        for value in calculated
    ]

    require(
        reportable_three
        == list(recorded["reportable_3"]),
        (
            f"{name} reportable interval differs: "
            f"{reportable_three} versus "
            f"{recorded['reportable_3']}"
        ),
    )

    row = summary_row(test_name)

    summary_interval = [
        float(value)
        for value in ast.literal_eval(
            str(row["null_mean"])
        )
    ]

    require(
        summary_interval == reportable_three,
        (
            f"{test_name} summary interval differs: "
            f"{summary_interval} versus {reportable_three}"
        ),
    )

    expected_observed = (
        observed["mechanism_consistency"]
        if metric == "mc"
        else observed["structural_agreement"]
    )

    require(
        close(
            float(row["observed"]),
            round(float(expected_observed), 4),
            tolerance=5e-5,
        ),
        f"{test_name} observed value differs",
    )


# ---------------------------------------------------------------------
# Ledger synchronization
# ---------------------------------------------------------------------

require(
    ledger.get("SA_total") == [99, 133],
    f"ledger SA_total is {ledger.get('SA_total')!r}",
)

require(
    close(
        ledger.get("SA"),
        99 / 133,
        tolerance=5e-5,
    ),
    f"ledger scalar SA is {ledger.get('SA')!r}",
)

require(
    ledger.get("SA_by_axis", {}).get("monomer", [])[:3]
    == [21, 28, 0.75],
    "ledger monomer-axis record differs from 21/28",
)

require(
    ledger.get("SA_by_axis", {}).get("fold", [])[:3]
    == [20, 26, 0.7692],
    "ledger complex-context record differs from 20/26",
)

require(
    ledger.get("SA_by_axis", {}).get("binding", [])[:3]
    == [24, 32, 0.75],
    "ledger binding record differs from 24/32",
)

require(
    ledger.get("SA_by_axis", {}).get("tier", [])[:3]
    == [34, 47, 0.7234],
    "ledger tier record differs from 34/47",
)

expected_cluster_sa = list(
    statistics["system_cluster_ci95"]["sa"]["rounded_4"]
)

require(
    ledger.get("SA_cluster_ci95")
    == expected_cluster_sa,
    (
        "ledger SA_cluster_ci95 differs from the verified draws: "
        f"{ledger.get('SA_cluster_ci95')} versus "
        f"{expected_cluster_sa}"
    ),
)

if "canonical_sha256" in ledger:
    require(
        ledger["canonical_sha256"]
        == expected_hashes["canonical_sha256"],
        "ledger canonical SHA-256 differs from the current canonical table",
    )


print("Canonical rows: 61")
print("Mechanism population: 57")
print("Primary structural agreement: 99/132")
print("All-row structural agreement: 99/133")
print(
    "Variant-bootstrap SA 95% CI:",
    statistics["variant_bootstrap_ci95"]["sa"]["reportable_3"],
)
print(
    "System-cluster SA 95% CI:",
    statistics["system_cluster_ci95"]["sa"]["reportable_3"],
)
print("STRESS-STATISTICS VERIFICATION: PASS")
