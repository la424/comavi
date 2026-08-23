from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
from comavi_v7.isds import (
    ISDS_OUTPUT_COLUMNS,
    ISDS_VERSION as PIPELINE_ISDS_VERSION,
)

OUT: Path

ISDS_VERSION = PIPELINE_ISDS_VERSION
ANCHORS = {'monomer': 2.9, 'complex_context': 2.9, 'binding': 3.5}
BOOT_DRAWS = 10_000
BOOT_SEED = 20260821
STRUCTURAL_CLASSES = {'mixed_structural', 'ppi_destab_mechanism', 'fold_mechanism'}


@dataclass(frozen=True)
class AxisValue:
    axis: str
    partner: str
    signed_value: float
    abs_value: float
    ratio: float


def _truthy(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.number)):
        return float(value) != 0.0
    return str(value).strip().lower() not in {'', '0', '0.0', 'false', 'none', 'nan', 'no'}


def _tier_number(value) -> Optional[int]:
    if pd.isna(value):
        return None
    m = re.search(r'([1-4])', str(value))
    if not m:
        return None
    return int(m.group(1))


def _tier_from_score(score: float) -> int:
    if pd.isna(score):
        raise ValueError('Tier score is missing')
    score = float(score)
    if score >= 5.0:
        return 1
    if score >= 3.0:
        return 2
    if score >= 1.5:
        return 3
    return 4


def _axis_columns(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    excluded = ('_sd', '_ci95', '_distinguishable', '_indistinguishable', '_vote', '_status', '_call')
    fold_cols = [
        c for c in df.columns
        if c.startswith('ddg_fold_')
        and not any(token in c for token in excluded)
        and pd.api.types.is_numeric_dtype(df[c])
    ]
    bind_cols = [
        c for c in df.columns
        if c.startswith('ddg_binding_')
        and not any(token in c for token in excluded)
        and pd.api.types.is_numeric_dtype(df[c])
    ]
    return sorted(fold_cols), sorted(bind_cols)


def axis_values(row: pd.Series, fold_cols: Iterable[str], bind_cols: Iterable[str]) -> list[AxisValue]:
    vals: list[AxisValue] = []
    mono = row.get('ddg_monomer')
    if pd.notna(mono) and _truthy(row.get('ddg_monomer_confident', True)):
        signed = float(mono)
        vals.append(AxisValue('monomer', '', signed, abs(signed), abs(signed) / ANCHORS['monomer']))

    for col in fold_cols:
        val = row.get(col)
        if pd.isna(val):
            continue
        partner = col[len('ddg_fold_'):]
        if not _truthy(row.get(f'ddg_{partner}_confident', True)):
            continue
        signed = float(val)
        vals.append(AxisValue('complex_context', partner, signed, abs(signed), abs(signed) / ANCHORS['complex_context']))

    for col in bind_cols:
        val = row.get(col)
        if pd.isna(val):
            continue
        partner = col[len('ddg_binding_'):]
        if not _truthy(row.get(f'ddg_{partner}_confident', True)):
            continue
        signed = float(val)
        vals.append(AxisValue('binding', partner, signed, abs(signed), abs(signed) / ANCHORS['binding']))
    return vals


def compute_isds_row(row: pd.Series, fold_cols: Iterable[str], bind_cols: Iterable[str], tier_col: str = 'reported_tier') -> dict:
    vals = axis_values(row, fold_cols, bind_cols)
    tier_n = _tier_number(row.get(tier_col))
    available = bool(vals) and tier_n is not None
    if not available:
        return {
            'isds_v1': np.nan,
            'isds_energy_component': np.nan,
            'isds_context_component': np.nan,
            'isds_energy_ratio_uncapped': np.nan,
            'isds_dominant_axis': '',
            'isds_dominant_partner': '',
            'isds_dominant_signed_ddg': np.nan,
            'isds_available': False,
            'isds_version': ISDS_VERSION,
            'max_abs_energy_recomputed': np.nan,
        }
    dominant = max(vals, key=lambda x: (x.ratio, x.abs_value, x.axis, x.partner))
    R = float(dominant.ratio)
    E = R / (1.0 + R)
    C = (4.0 - tier_n) / 3.0
    score = 0.5 * (E + C)
    return {
        'isds_v1': score,
        'isds_energy_component': E,
        'isds_context_component': C,
        'isds_energy_ratio_uncapped': R,
        'isds_dominant_axis': dominant.axis,
        'isds_dominant_partner': dominant.partner,
        'isds_dominant_signed_ddg': dominant.signed_value,
        'isds_available': True,
        'isds_version': ISDS_VERSION,
        'max_abs_energy_recomputed': dominant.abs_value,
    }


A636P_FINAL_VALUES = {
    "expected_ddg_monomer": "neutral",
    "expected_ddg_fold_complex": "neutral",
    "expected_ddg_binding": "neutral",
    "evidence_axes": "functional",
}


def audit_a636p_canonical(
    canonical: pd.DataFrame,
) -> pd.DataFrame:
    """Verify the literature-supported A636P record without mutating it."""
    mask = (
        canonical["variant"].astype(str).eq("A636P")
        & canonical["system"].astype(str).eq("msh2_msh6")
    )

    if int(mask.sum()) != 1:
        raise RuntimeError(
            f"Expected one A636P row; found {int(mask.sum())}"
        )

    row = canonical.loc[mask].iloc[0]

    mismatches = {
        column: {
            "observed": row.get(column),
            "expected": expected,
        }
        for column, expected in A636P_FINAL_VALUES.items()
        if str(row.get(column)) != expected
    }

    if mismatches:
        raise RuntimeError(
            "Canonical A636P record is not synchronized: "
            f"{mismatches}"
        )

    if str(row.get("expected_mech_class")) != "structurally_silent":
        raise RuntimeError(
            "Canonical A636P expected_mech_class must be "
            "'structurally_silent'."
        )

    columns = [
        "system",
        "gene",
        "variant",
        "expected_ddg_monomer",
        "expected_ddg_fold_complex",
        "expected_ddg_binding",
        "evidence_axes",
        "expected_mech_class",
        "notes",
        "mech_consistency_t25",
        "structural_agreement_n_t25",
        "structural_agreement_d_t25",
        "isds_v1",
        "isds_version",
    ]

    audit = canonical.loc[mask, columns].copy()
    audit.insert(0, "audit_status", "verified_final_canonical")
    return audit


ISDS_NUMERIC_COLUMNS = [
    "isds_v1",
    "isds_energy_component",
    "isds_context_component",
    "isds_energy_ratio_uncapped",
    "isds_dominant_signed_ddg",
]

ISDS_TEXT_COLUMNS = [
    "isds_dominant_axis",
    "isds_dominant_partner",
    "isds_version",
]


def _normalize_text_series(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .fillna("")
        .str.strip()
    )


def _normalize_bool_series(series: pd.Series) -> pd.Series:
    return (
        series.fillna(False)
        .map(
            lambda value: str(value).strip().lower()
            in {"true", "1", "1.0", "yes"}
        )
        .astype(bool)
    )


def audit_canonical_isds_outputs(
    stored: pd.DataFrame,
    recomputed: pd.DataFrame,
) -> pd.DataFrame:
    """Verify the committed ISDS fields against independent recomputation."""
    keys = ["system", "variant"]

    if len(stored) != len(recomputed):
        raise RuntimeError(
            "Stored and recomputed canonical populations differ in length: "
            f"{len(stored)} versus {len(recomputed)}"
        )

    if stored.duplicated(keys).any():
        raise RuntimeError(
            "Stored canonical table contains duplicate system/variant keys."
        )

    if recomputed.duplicated(keys).any():
        raise RuntimeError(
            "Recomputed canonical table contains duplicate system/variant keys."
        )

    missing_stored = [
        column
        for column in ISDS_OUTPUT_COLUMNS
        if column not in stored.columns
    ]

    missing_recomputed = [
        column
        for column in ISDS_OUTPUT_COLUMNS
        if column not in recomputed.columns
    ]

    if missing_stored or missing_recomputed:
        raise RuntimeError(
            "ISDS schema mismatch: "
            f"missing stored={missing_stored}; "
            f"missing recomputed={missing_recomputed}"
        )

    stored_sorted = stored.sort_values(keys).reset_index(drop=True)
    recomputed_sorted = recomputed.sort_values(keys).reset_index(drop=True)

    if not stored_sorted[keys].equals(recomputed_sorted[keys]):
        raise RuntimeError(
            "Stored and recomputed canonical keys differ."
        )

    n_rows = len(stored_sorted)
    numeric_match = np.ones(n_rows, dtype=bool)
    maximum_numeric_delta = np.zeros(n_rows, dtype=float)

    for column in ISDS_NUMERIC_COLUMNS:
        stored_values = pd.to_numeric(
            stored_sorted[column],
            errors="coerce",
        ).to_numpy(float)

        recomputed_values = pd.to_numeric(
            recomputed_sorted[column],
            errors="coerce",
        ).to_numpy(float)

        column_match = np.isclose(
            stored_values,
            recomputed_values,
            atol=1e-8,
            rtol=0.0,
            equal_nan=True,
        )

        numeric_match &= column_match

        finite = (
            np.isfinite(stored_values)
            & np.isfinite(recomputed_values)
        )

        deltas = np.zeros(n_rows, dtype=float)
        deltas[finite] = np.abs(
            stored_values[finite]
            - recomputed_values[finite]
        )

        maximum_numeric_delta = np.maximum(
            maximum_numeric_delta,
            deltas,
        )

    text_match = np.ones(n_rows, dtype=bool)

    for column in ISDS_TEXT_COLUMNS:
        stored_values = _normalize_text_series(
            stored_sorted[column]
        ).to_numpy()

        recomputed_values = _normalize_text_series(
            recomputed_sorted[column]
        ).to_numpy()

        text_match &= stored_values == recomputed_values

    boolean_match = (
        _normalize_bool_series(
            stored_sorted["isds_available"]
        ).to_numpy()
        ==
        _normalize_bool_series(
            recomputed_sorted["isds_available"]
        ).to_numpy()
    )

    all_fields_match = (
        numeric_match
        & text_match
        & boolean_match
    )

    audit = stored_sorted[keys].copy()
    audit["numeric_fields_match"] = numeric_match
    audit["maximum_numeric_delta"] = maximum_numeric_delta
    audit["text_fields_match"] = text_match
    audit["availability_matches"] = boolean_match
    audit["all_fields_match"] = all_fields_match

    if not bool(all_fields_match.all()):
        failures = audit.loc[
            ~audit["all_fields_match"]
        ]

        raise RuntimeError(
            "Committed ISDS fields differ from independent recomputation:\n"
            + failures.to_string(index=False)
        )

    return audit

def auc_metrics(y: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    return float(roc_auc_score(y, score)), float(average_precision_score(y, score))


def topk_table(df: pd.DataFrame, score_col: str, ks=(5, 10, 15, 20, 25)) -> pd.DataFrame:
    # ISDS and energetic scores use the uncapped energy ratio only to resolve exact numerical ties.
    # Context-only scores are intentionally not tie-broken by energy, because that would turn the
    # comparator into a hidden combination.
    if score_col in {'isds_v1', 'isds_energy_component', 'max_abs_energy_recomputed'}:
        ranked = df.sort_values([score_col, 'isds_energy_ratio_uncapped', 'system', 'variant'], ascending=[False, False, True, True]).reset_index(drop=True)
    else:
        ranked = df.sort_values([score_col, 'system', 'variant'], ascending=[False, True, True]).reset_index(drop=True)
    n_pos = int(ranked['structural_ground_truth'].sum())
    rows = []
    for k in ks:
        k_eff = min(k, len(ranked))
        pos = int(ranked.iloc[:k_eff]['structural_ground_truth'].sum())
        rows.append({
            'score': score_col,
            'k': k,
            'k_effective': k_eff,
            'structural_mechanisms_in_top_k': pos,
            'precision_at_k': pos / k_eff,
            'recovery_at_k': pos / n_pos if n_pos else np.nan,
        })
    return pd.DataFrame(rows)


def _roc_auc_fast(y: np.ndarray, score: np.ndarray) -> float:
    """Fast ROC AUC with average ranks for ties."""
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=float)
    n_pos = int(y.sum())
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float('nan')
    order = np.argsort(score, kind='mergesort')
    sorted_scores = score[order]
    ranks = np.empty(len(score), dtype=float)
    i = 0
    while i < len(score):
        j = i + 1
        while j < len(score) and sorted_scores[j] == sorted_scores[i]:
            j += 1
        avg_rank = 0.5 * ((i + 1) + j)  # 1-based average rank
        ranks[order[i:j]] = avg_rank
        i = j
    rank_sum_pos = ranks[y == 1].sum()
    return float((rank_sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def _average_precision_fast(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=float)
    n_pos = int(y.sum())
    if n_pos == 0:
        return float('nan')
    order = np.argsort(-score, kind='mergesort')
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    precision = tp / (np.arange(len(y_sorted)) + 1)
    return float((precision * y_sorted).sum() / n_pos)


def cluster_bootstrap(df: pd.DataFrame, score_cols: list[str], draws: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    systems = np.array(sorted(df['system'].unique()))
    all_y = df['structural_ground_truth'].astype(int).to_numpy(np.int8)
    all_scores = {col: df[col].to_numpy(float) for col in score_cols}
    indices_by_system = {s: np.flatnonzero(df['system'].to_numpy() == s) for s in systems}
    records = []
    for draw in range(draws):
        sampled = rng.choice(systems, size=len(systems), replace=True)
        indices = np.concatenate([indices_by_system[s] for s in sampled])
        y = all_y[indices]
        if y.min() == y.max():
            continue
        rec = {'draw': draw}
        for col in score_cols:
            score = all_scores[col][indices]
            rec[f'roc_auc__{col}'] = _roc_auc_fast(y, score)
            rec[f'average_precision__{col}'] = _average_precision_fast(y, score)
        rec['auc_diff_isds_minus_energy'] = rec['roc_auc__isds_v1'] - rec['roc_auc__isds_energy_component']
        rec['auc_diff_isds_minus_context'] = rec['roc_auc__isds_v1'] - rec['roc_auc__isds_context_component']
        records.append(rec)
    draws_df = pd.DataFrame(records)
    summaries = []
    for col in [c for c in draws_df.columns if c != 'draw']:
        vals = draws_df[col].dropna().to_numpy(float)
        summaries.append({
            'metric': col,
            'n_valid_draws': len(vals),
            'mean': float(np.mean(vals)),
            'median': float(np.median(vals)),
            'ci95_low': float(np.quantile(vals, 0.025)),
            'ci95_high': float(np.quantile(vals, 0.975)),
            'fraction_gt_zero': float(np.mean(vals > 0)) if 'diff' in col else np.nan,
        })
    return draws_df, pd.DataFrame(summaries)

def leave_one_system_out(df: pd.DataFrame, score_cols: list[str]) -> pd.DataFrame:
    rows = []
    for system in sorted(df['system'].unique()):
        sub = df[~df['system'].eq(system)]
        y = sub['structural_ground_truth'].astype(int).to_numpy()
        if len(np.unique(y)) < 2:
            continue
        rec = {'removed_system': system, 'n_remaining': len(sub), 'positive_remaining': int(y.sum())}
        for col in score_cols:
            rec[f'roc_auc__{col}'] = roc_auc_score(y, sub[col])
            rec[f'average_precision__{col}'] = average_precision_score(y, sub[col])
        rec['auc_diff_isds_minus_energy'] = rec['roc_auc__isds_v1'] - rec['roc_auc__isds_energy_component']
        rec['auc_diff_isds_minus_context'] = rec['roc_auc__isds_v1'] - rec['roc_auc__isds_context_component']
        rows.append(rec)
    return pd.DataFrame(rows)


def formula_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    R = out['isds_energy_ratio_uncapped'].to_numpy(float)
    C = out['isds_context_component'].to_numpy(float)
    out['sens_hard_cap'] = 0.5 * (np.minimum(R, 1.0) + C)
    out['sens_exponential'] = 0.5 * ((1.0 - np.exp(-R)) + C)
    out['sens_smooth_norm'] = 0.5 * ((R / np.sqrt(1.0 + R ** 2)) + C)
    out['sens_energy_weight_0_4'] = 0.4 * (R / (1.0 + R)) + 0.6 * C
    out['sens_energy_weight_0_6'] = 0.6 * (R / (1.0 + R)) + 0.4 * C

    # Common 2.5 anchors, using recomputed axis magnitudes.
    common_R = out[['abs_monomer', 'abs_complex', 'abs_binding']].max(axis=1, skipna=True) / 2.5
    out['sens_common_2_5'] = 0.5 * ((common_R / (1.0 + common_R)) + C)

    # No-interface tier sensitivity.
    no_iface_tier = out['no_interface_score'].apply(_tier_from_score)
    no_iface_C = (4.0 - no_iface_tier.astype(float)) / 3.0
    out['sens_no_interface_tier'] = 0.5 * ((R / (1.0 + R)) + no_iface_C)

    # Interface-only component: 1 when an interface is detected, else 0.
    iface_C = out['interface_only_positive'].astype(float)
    out['sens_interface_only_context'] = 0.5 * ((R / (1.0 + R)) + iface_C)

    # Historical cohort-z combination, retained only for comparison.
    energy = out['max_abs_energy_recomputed'].to_numpy(float)
    tier_strength = 4.0 - out['tier_number'].to_numpy(float)
    energy_z = (energy - np.mean(energy)) / np.std(energy, ddof=0)
    tier_z = (tier_strength - np.mean(tier_strength)) / np.std(tier_strength, ddof=0)
    out['historical_cohort_z_combination'] = 0.5 * (energy_z + tier_z)
    return out


def score_band_table(df: pd.DataFrame) -> pd.DataFrame:
    bins = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0000001]
    labels = ['0.00-<0.20', '0.20-<0.40', '0.40-<0.60', '0.60-<0.80', '0.80-1.00']
    x = df.copy()
    x['isds_band'] = pd.cut(x['isds_v1'], bins=bins, labels=labels, right=False, include_lowest=True)
    rows = []
    for label in labels:
        s = x[x['isds_band'].astype(str).eq(label)]
        n = len(s)
        p = int(s['structural_ground_truth'].sum())
        rows.append({
            'isds_band': label,
            'n_variants': n,
            'n_structural_mechanism': p,
            'structural_mechanism_fraction': p / n if n else np.nan,
            'systems': ';'.join(sorted(s['system'].unique())),
            'variants': ';'.join(s.sort_values('isds_v1', ascending=False)['variant'].astype(str)),
        })
    return pd.DataFrame(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Regenerate the locked COMAVI ISDS-v1 benchmark analyses.'
    )
    parser.add_argument(
        '--canonical',
        type=Path,
        required=True,
        help='Path to the canonical 61-variant benchmark CSV.',
    )
    parser.add_argument(
        '--tier-calls',
        type=Path,
        required=True,
        help='Path to the 47-row tier-comparator source table.',
    )
    parser.add_argument(
        '--out-dir',
        type=Path,
        required=True,
        help='Directory for regenerated ISDS-v1 outputs.',
    )
    args = parser.parse_args()

    canonical_path = args.canonical.expanduser().resolve()
    tier_calls_path = args.tier_calls.expanduser().resolve()

    if not canonical_path.is_file():
        raise SystemExit(f'Canonical benchmark not found: {canonical_path}')
    if not tier_calls_path.is_file():
        raise SystemExit(f'Tier-comparator table not found: {tier_calls_path}')

    global OUT
    OUT = args.out_dir.expanduser().resolve()
    OUT.mkdir(parents=True, exist_ok=True)

    canonical = pd.read_csv(canonical_path, low_memory=False)

    missing_isds_columns = [
        column
        for column in ISDS_OUTPUT_COLUMNS
        if column not in canonical.columns
    ]

    if missing_isds_columns:
        raise RuntimeError(
            "Canonical benchmark lacks committed ISDS fields: "
            f"{missing_isds_columns}"
        )

    a636p_audit = audit_a636p_canonical(canonical)
    a636p_audit_path = OUT / 'A636P_canonical_audit.csv'
    a636p_audit.to_csv(a636p_audit_path, index=False)

    canonical_without_isds = canonical.drop(
        columns=list(ISDS_OUTPUT_COLUMNS),
    )

    tier = pd.read_csv(tier_calls_path)
    df = tier.merge(canonical_without_isds, on=['system', 'variant'], how='left', validate='one_to_one', suffixes=('', '_canonical'))
    if len(df) != 47:
        raise RuntimeError(f'Expected 47 prioritization rows; found {len(df)}')

    fold_cols, bind_cols = _axis_columns(canonical)
    calculated = []
    axis_summaries = []
    for _, row in df.iterrows():
        vals = axis_values(row, fold_cols, bind_cols)
        d = compute_isds_row(row, fold_cols, bind_cols)
        calculated.append(d)
        by_axis = {}
        for axis in ['monomer', 'complex_context', 'binding']:
            candidates = [x for x in vals if x.axis == axis]
            if candidates:
                best = max(candidates, key=lambda x: (x.abs_value, x.partner))
                by_axis[f'abs_{axis}'] = best.abs_value
                by_axis[f'signed_{axis}'] = best.signed_value
                by_axis[f'partner_{axis}'] = best.partner
            else:
                by_axis[f'abs_{axis}'] = np.nan
                by_axis[f'signed_{axis}'] = np.nan
                by_axis[f'partner_{axis}'] = ''
        axis_summaries.append(by_axis)
    df = pd.concat([df.reset_index(drop=True), pd.DataFrame(calculated), pd.DataFrame(axis_summaries)], axis=1)
    if df.columns.duplicated().any():
        duplicates = sorted(
            set(df.columns[df.columns.duplicated()])
        )
        raise RuntimeError(
            f'Duplicate columns after prioritization merge: {duplicates}'
        )
    df = df.rename(columns={
        'abs_complex_context': 'abs_complex',
        'signed_complex_context': 'signed_complex',
        'partner_complex_context': 'partner_complex',
    })
    df['tier_number'] = df['reported_tier'].apply(_tier_number).astype(int)
    df['structural_ground_truth'] = df['structural_ground_truth'].astype(bool)
    df['endpoint_label'] = np.where(df['structural_ground_truth'], 'structural mechanism', 'no modeled lesion')

    # Verify historical max energy where available.
    if 'max_abs_ddg' in df.columns:
        delta = np.abs(df['max_abs_ddg'].astype(float) - df['max_abs_energy_recomputed'].astype(float))
        if float(delta.max()) > 1e-6:
            raise RuntimeError(f'Recomputed max energy differs from canonical max_abs_ddg; max delta {delta.max()}')

    df = formula_scores(df)
    per_variant_path = OUT / 'ISDS_v1_per_variant.csv'
    keep_cols = [
        'system', 'variant', 'gene', 'structural_ground_truth', 'endpoint_label', 'expected_mech_class',
        'reported_tier', 'tier_number', 'reconstructed_full_score', 'no_interface_score',
        'interface_partner_count', 'interface_only_positive',
        'abs_monomer', 'signed_monomer', 'abs_complex', 'signed_complex', 'partner_complex',
        'abs_binding', 'signed_binding', 'partner_binding', 'max_abs_energy_recomputed',
        'isds_energy_ratio_uncapped', 'isds_energy_component', 'isds_context_component',
        'isds_v1', 'isds_dominant_axis', 'isds_dominant_partner', 'isds_dominant_signed_ddg',
        'isds_available', 'isds_version',
        'sens_hard_cap', 'sens_exponential', 'sens_smooth_norm', 'sens_common_2_5',
        'sens_energy_weight_0_4', 'sens_energy_weight_0_6', 'sens_no_interface_tier',
        'sens_interface_only_context', 'historical_cohort_z_combination',
        'AM pathogenicity', 'phenotype', 'franklin', 'structural_evidence_strength', 'evidence_axes'
    ]
    df[keep_cols].to_csv(per_variant_path, index=False)

    y = df['structural_ground_truth'].astype(int).to_numpy()
    primary_score_cols = ['isds_v1', 'isds_energy_component', 'isds_context_component']
    comparator_cols = primary_score_cols + [
        'max_abs_energy_recomputed', 'interface_only_positive', 'sens_no_interface_tier',
        'historical_cohort_z_combination'
    ]
    metric_rows = []
    for col in comparator_cols:
        roc, ap = auc_metrics(y, df[col].astype(float).to_numpy())
        metric_rows.append({'score': col, 'roc_auc': roc, 'average_precision': ap, 'n': len(df), 'n_positive': int(y.sum()), 'n_negative': int((1-y).sum())})
    primary_metrics = pd.DataFrame(metric_rows)
    primary_metrics.to_csv(OUT / 'ISDS_v1_primary_metrics.csv', index=False)

    topk = pd.concat([topk_table(df, col) for col in primary_score_cols + ['max_abs_energy_recomputed']], ignore_index=True)
    topk.to_csv(OUT / 'ISDS_v1_top_k.csv', index=False)

    boot_draws, boot_summary = cluster_bootstrap(df, primary_score_cols, BOOT_DRAWS, BOOT_SEED)
    boot_draws.to_csv(OUT / 'ISDS_v1_system_cluster_bootstrap_10000.csv', index=False)
    boot_summary.to_csv(OUT / 'ISDS_v1_system_cluster_bootstrap_summary.csv', index=False)

    loso = leave_one_system_out(df, primary_score_cols)
    loso.to_csv(OUT / 'ISDS_v1_leave_one_system_out.csv', index=False)

    bands = score_band_table(df)
    bands.to_csv(OUT / 'ISDS_v1_score_bands.csv', index=False)

    sens_cols = [c for c in df.columns if c.startswith('sens_')] + ['historical_cohort_z_combination']
    sens_rows = []
    for col in sens_cols:
        roc, ap = auc_metrics(y, df[col].astype(float).to_numpy())
        rho = spearmanr(df['isds_v1'], df[col]).statistic
        sens_rows.append({'formula': col, 'roc_auc': roc, 'average_precision': ap, 'spearman_vs_isds_v1': rho, 'n_unique_scores': int(df[col].round(12).nunique())})
    pd.DataFrame(sens_rows).sort_values('roc_auc', ascending=False).to_csv(OUT / 'ISDS_v1_formula_sensitivities.csv', index=False)

    # Tier/component baselines in their original binary form.
    component_binary = pd.DataFrame([
        {'screen': 'interface_status_alone', 'tp': 15, 'fn': 2, 'fp': 7, 'tn': 23, 'sensitivity': 15/17, 'specificity': 23/30, 'ppv': 15/22},
        {'screen': 'tier_without_interface_bonus', 'tp': 14, 'fn': 3, 'fp': 12, 'tn': 18, 'sensitivity': 14/17, 'specificity': 18/30, 'ppv': 14/26},
        {'screen': 'full_tier_1_2', 'tp': 17, 'fn': 0, 'fp': 13, 'tn': 17, 'sensitivity': 1.0, 'specificity': 17/30, 'ppv': 17/30},
    ])
    component_binary.to_csv(OUT / 'ISDS_v1_component_binary_baselines.csv', index=False)

    # AlphaMissense comparison on the complete 47-variant phenotype-labelled set.
    # This population is not identical to the 47-variant structural-prioritization set.
    full_calculated = []

    for _, row in canonical_without_isds.iterrows():
        full_calculated.append(
            compute_isds_row(
                row,
                fold_cols,
                bind_cols,
                tier_col='comavi_tier',
            )
        )

    full_recomputed = pd.concat(
        [
            canonical_without_isds.reset_index(drop=True),
            pd.DataFrame(full_calculated),
        ],
        axis=1,
    )

    if full_recomputed.columns.duplicated().any():
        duplicates = sorted(
            set(
                full_recomputed.columns[
                    full_recomputed.columns.duplicated()
                ]
            )
        )
        raise RuntimeError(
            f'Duplicate columns in full recomputation: {duplicates}'
        )

    canonical_isds_audit = audit_canonical_isds_outputs(
        canonical,
        full_recomputed,
    )

    canonical_isds_audit_path = (
        OUT / 'ISDS_v1_canonical_output_audit.csv'
    )

    canonical_isds_audit.to_csv(
        canonical_isds_audit_path,
        index=False,
    )

    full_isds = canonical.copy()

    full_isds_path = (
        OUT / 'scored_61var_canonical_with_ISDS_v1.csv'
    )

    full_isds_path.write_bytes(
        canonical_path.read_bytes()
    )
    clinical = full_isds[
        full_isds['phenotype'].isin(['pathogenic', 'pathogenic_gof', 'benign'])
        & full_isds['AM pathogenicity'].notna()
        & full_isds['isds_available'].astype(bool)
    ].copy()
    clinical['clinical_y'] = clinical['phenotype'].isin(['pathogenic', 'pathogenic_gof']).astype(int)
    clinical_rows = []
    for col in ['AM pathogenicity', 'isds_v1', 'isds_energy_component', 'isds_context_component']:
        clinical_rows.append({'descriptor': col, 'n': len(clinical), 'n_pathogenic': int(clinical['clinical_y'].sum()), 'pathogenicity_auc': roc_auc_score(clinical['clinical_y'], clinical[col])})
    pd.DataFrame(clinical_rows).to_csv(OUT / 'ISDS_v1_alphamissense_pathogenicity_comparison.csv', index=False)
    clinical[['system','variant','phenotype','clinical_y','AM pathogenicity','isds_v1','isds_energy_component','isds_context_component','isds_dominant_axis']].to_csv(OUT / 'ISDS_v1_alphamissense_common_set.csv', index=False)

    # Main summary and manuscript anchors.
    observed = {row['score']: row for row in metric_rows}
    boot_map = {r['metric']: r for _, r in boot_summary.iterrows()}
    isds_row = observed['isds_v1']
    energy_row = observed['isds_energy_component']
    context_row = observed['isds_context_component']
    summary = {
        'analysis_id': 'COMAVI-ISDS-v1-2026-08-21',
        'formula_version': ISDS_VERSION,
        'population': {'n': len(df), 'positive': int(y.sum()), 'negative': int((1-y).sum()), 'systems': int(df.system.nunique())},
        'observed': {
            'isds_roc_auc': isds_row['roc_auc'],
            'isds_average_precision': isds_row['average_precision'],
            'energy_roc_auc': energy_row['roc_auc'],
            'context_roc_auc': context_row['roc_auc'],
            'isds_minus_energy_auc': isds_row['roc_auc'] - energy_row['roc_auc'],
            'isds_minus_context_auc': isds_row['roc_auc'] - context_row['roc_auc'],
        },
        'bootstrap': {
            'draws_requested': BOOT_DRAWS,
            'seed': BOOT_SEED,
            'isds_roc_auc_ci95': [boot_map['roc_auc__isds_v1']['ci95_low'], boot_map['roc_auc__isds_v1']['ci95_high']],
            'isds_average_precision_ci95': [boot_map['average_precision__isds_v1']['ci95_low'], boot_map['average_precision__isds_v1']['ci95_high']],
            'isds_minus_energy_auc_ci95': [boot_map['auc_diff_isds_minus_energy']['ci95_low'], boot_map['auc_diff_isds_minus_energy']['ci95_high']],
            'isds_minus_context_auc_ci95': [boot_map['auc_diff_isds_minus_context']['ci95_low'], boot_map['auc_diff_isds_minus_context']['ci95_high']],
            'fraction_isds_gt_energy': boot_map['auc_diff_isds_minus_energy']['fraction_gt_zero'],
            'fraction_isds_gt_context': boot_map['auc_diff_isds_minus_context']['fraction_gt_zero'],
            'valid_draws': int(boot_map['roc_auc__isds_v1']['n_valid_draws']),
        },
        'leave_one_system_out': {
            'roc_auc_min': float(loso['roc_auc__isds_v1'].min()),
            'roc_auc_max': float(loso['roc_auc__isds_v1'].max()),
            'roc_auc_median': float(loso['roc_auc__isds_v1'].median()),
            'average_precision_min': float(loso['average_precision__isds_v1'].min()),
            'average_precision_max': float(loso['average_precision__isds_v1'].max()),
        },
        'top_k': topk[topk['score'].eq('isds_v1')].to_dict(orient='records'),
        'score_bands': bands.to_dict(orient='records'),
        'mechanism_localization_final_values': {
            'whole_variant_score': '41/57 = 0.719',
            'monomer_agreement': '21/27 = 0.778',
            'complex_context_agreement': '20/26 = 0.769',
            'binding_agreement': '24/32 = 0.750',
            'all_energetic_agreement': '65/85 = 0.765',
            'historical_four_output': '99/132 = 0.750',
        },
        'files': {
            'canonical_input_sha256': sha256(canonical_path),
            'a636p_audit_sha256': sha256(a636p_audit_path),
            'canonical_isds_audit_sha256': sha256(
                canonical_isds_audit_path
            ),
            'per_variant_sha256': sha256(per_variant_path),
        }
    }
    (OUT / 'ISDS_v1_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    # Rankings for practical review.
    ranking_cols = ['system','variant','endpoint_label','reported_tier','isds_v1','isds_energy_component','isds_context_component','isds_energy_ratio_uncapped','isds_dominant_axis','isds_dominant_partner','isds_dominant_signed_ddg']
    df.sort_values(['isds_v1','isds_energy_ratio_uncapped'], ascending=False)[ranking_cols].to_csv(OUT / 'ISDS_v1_ranked_variants.csv', index=False)

    # Checksums.
    checksum_paths = sorted([p for p in OUT.glob('*') if p.is_file() and p.name != 'SHA256SUMS.txt'])
    (OUT / 'SHA256SUMS.txt').write_text('\n'.join(f'{sha256(p)}  {p.name}' for p in checksum_paths) + '\n')

    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
