"""Cohort-independent Integrated Structural-Disruption Score (ISDS-v1).

ISDS-v1 combines the strongest axis-normalized absolute FoldX magnitude with
COMAVI's ordinal structural-context tier.  It is a prioritization index, not a
pathogenicity probability and not a replacement for the signed mechanism
profile.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Mapping, Optional, Sequence

import pandas as pd

ISDS_VERSION = "ISDS-v1"
ISDS_OUTPUT_COLUMNS = (
    "isds_v1",
    "isds_energy_component",
    "isds_context_component",
    "isds_energy_ratio_uncapped",
    "isds_dominant_axis",
    "isds_dominant_partner",
    "isds_dominant_signed_ddg",
    "isds_available",
    "isds_version",
)
ISDS_ENERGY_ANCHORS = {
    "monomer": 2.9,
    "complex_context": 2.9,
    "binding": 3.5,
}


@dataclass(frozen=True)
class ISDSAxisEvidence:
    axis: str
    partner: str
    signed_ddg: float
    ratio: float


def _truthy(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value) != 0.0
    return str(value).strip().lower() not in {"", "0", "0.0", "false", "none", "nan", "no"}


def _tier_number(value) -> Optional[int]:
    if pd.isna(value):
        return None
    match = re.search(r"([1-4])", str(value))
    return int(match.group(1)) if match else None


def _valid_number(value) -> bool:
    return value is not None and not pd.isna(value)


def _collect_axis_evidence(row: pd.Series, partner_labels: Sequence[str]) -> list[ISDSAxisEvidence]:
    evidence: list[ISDSAxisEvidence] = []

    monomer = row.get("ddg_monomer")
    if _valid_number(monomer) and _truthy(row.get("ddg_monomer_confident", True)):
        signed = float(monomer)
        evidence.append(ISDSAxisEvidence(
            axis="monomer",
            partner="",
            signed_ddg=signed,
            ratio=abs(signed) / ISDS_ENERGY_ANCHORS["monomer"],
        ))

    for partner in partner_labels:
        confidence = row.get(f"ddg_{partner}_confident", True)
        if not _truthy(confidence):
            continue

        fold_value = row.get(f"ddg_fold_{partner}")
        if _valid_number(fold_value):
            signed = float(fold_value)
            evidence.append(ISDSAxisEvidence(
                axis="complex_context",
                partner=str(partner),
                signed_ddg=signed,
                ratio=abs(signed) / ISDS_ENERGY_ANCHORS["complex_context"],
            ))

        binding_value = row.get(f"ddg_binding_{partner}")
        if _valid_number(binding_value):
            signed = float(binding_value)
            evidence.append(ISDSAxisEvidence(
                axis="binding",
                partner=str(partner),
                signed_ddg=signed,
                ratio=abs(signed) / ISDS_ENERGY_ANCHORS["binding"],
            ))

    return evidence


def calculate_isds_v1(row: pd.Series, partner_labels: Sequence[str]) -> pd.Series:
    """Calculate ISDS-v1 and its transparent components for one variant row.

    Formula
    -------
    R = max(|DDG_monomer|/2.9, |DDG_complex|/2.9, |DDG_binding|/3.5)
    E = R/(1+R)
    C = (4-tier_number)/3
    ISDS-v1 = (E+C)/2

    Only valid, confidence-passing energetic axes enter R.  ISDS is unavailable
    when no energetic axis is valid or when the structural-context tier is
    unavailable.  The 1.0 kcal/mol binary mechanism threshold is intentionally
    not applied inside the continuous score.
    """
    tier_n = _tier_number(row.get("comavi_tier"))
    evidence = _collect_axis_evidence(row, partner_labels)
    available = tier_n is not None and bool(evidence)

    base = {
        "isds_v1": pd.NA,
        "isds_energy_component": pd.NA,
        "isds_context_component": pd.NA,
        "isds_energy_ratio_uncapped": pd.NA,
        "isds_dominant_axis": "",
        "isds_dominant_partner": "",
        "isds_dominant_signed_ddg": pd.NA,
        "isds_available": bool(available),
        "isds_version": ISDS_VERSION,
    }
    if not available:
        return pd.Series(base)

    dominant = max(evidence, key=lambda x: (x.ratio, abs(x.signed_ddg), x.axis, x.partner))
    ratio = float(dominant.ratio)
    energy_component = ratio / (1.0 + ratio)
    context_component = (4.0 - float(tier_n)) / 3.0
    score = 0.5 * (energy_component + context_component)

    base.update({
        "isds_v1": round(score, 8),
        "isds_energy_component": round(energy_component, 8),
        "isds_context_component": round(context_component, 8),
        "isds_energy_ratio_uncapped": round(ratio, 8),
        "isds_dominant_axis": dominant.axis,
        "isds_dominant_partner": dominant.partner,
        "isds_dominant_signed_ddg": round(float(dominant.signed_ddg), 8),
    })
    return pd.Series(base)


def add_isds_v1_columns(df: pd.DataFrame, partner_labels: Sequence[str]) -> pd.DataFrame:
    """Return a copy of *df* with versioned ISDS-v1 output columns appended."""
    out = df.copy()
    isds = out.apply(lambda row: calculate_isds_v1(row, partner_labels), axis=1)
    for column in isds.columns:
        out[column] = isds[column]
    return out


def infer_partner_labels(columns: Iterable[str]) -> list[str]:
    """Infer partner labels from existing ddg_fold_* and ddg_binding_* columns."""
    excluded_suffixes = (
        "_sd", "_ci95_internal_low", "_ci95_internal_high",
        "_ci95_sapozhnikov_low", "_ci95_sapozhnikov_high",
        "_distinguishable_internal_from_0", "_distinguishable_sapozhnikov_from_0",
        "_indistinguishable", "_vote_strict", "_vote_relaxed",
    )
    labels = set()
    for column in columns:
        prefix = None
        if column.startswith("ddg_fold_"):
            prefix = "ddg_fold_"
        elif column.startswith("ddg_binding_"):
            prefix = "ddg_binding_"
        if prefix is None or column.endswith(excluded_suffixes):
            continue
        labels.add(column[len(prefix):])
    return sorted(labels)
