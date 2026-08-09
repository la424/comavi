"""Does the structural-evidence tier add discriminative signal to FoldX energy?

The tier ranks pathogenicity (§3.9) and the energy axes detect mechanism
(§3.2-3.4). §3.11 asks whether the two interact. This script answers the
narrower and cleanly testable version of that question: as *pathogenicity
classifiers*, does combining them beat energy alone?

Both predictors are put on one scale (z-scored, tier inverted so that stronger
structural evidence scores higher) and summed with equal weight. Equal weighting
is deliberate: any fitted weight would be fit and evaluated on the same 49
variants, and the resulting AUC would be optimistic by an amount this benchmark
cannot estimate.

The uncertainty is a system-level cluster bootstrap, matching the paper's
primary interval convention -- variants within a complex are not independent.

Writes reference_outputs/COMAVI_tier_auc_decomposition.json.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
N_BOOT = 2000
SEED = 0


def _load_concordance():
    spec = importlib.util.spec_from_file_location("ac", REPO / "scripts/apply_concordance_v5.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ac"] = mod
    spec.loader.exec_module(mod)
    return mod


def auc(score, pos) -> float:
    """Mann-Whitney AUC. Ties get midranks, so this is exact under ties."""
    score, pos = np.asarray(score, float), np.asarray(pos, bool)
    r = stats.rankdata(score)
    n1, n0 = pos.sum(), (~pos).sum()
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _norm_label(s) -> str:
    s = str(s).lower()
    if "conflict" in s or "uncertain" in s:
        return "vus"
    if "pathogenic" in s:
        return "pathogenic"
    if "benign" in s:
        return "benign"
    return "other"


def load_labeled() -> pd.DataFrame:
    """Benchmark rows carrying both a clinical label and a tier.

    `phenotype` is null for the brca1_brct rows -- those variants' clinical
    calls live in the BRCT supplement's clinvar_germline column. Without this
    fallback the pathogenic count is silently wrong.
    """
    ac = _load_concordance()
    sc = pd.read_csv(REPO / "reference_outputs/scored_61var_canonical.csv")
    partners = ac.discover_partners(sc)

    brct = pd.read_csv(REPO / "supplement/brct/brct_foldx_concordance.csv")
    bmap = dict(zip(brct["variant"], brct["clinvar_germline"].astype(str)))

    def label_of(r):
        raw = r.get("phenotype")
        if (pd.isna(raw) or str(raw).strip() == "") and r["system"] == "brca1_brct":
            raw = bmap.get(r["variant"])
        return _norm_label(raw)

    sc["label"] = sc.apply(label_of, axis=1)
    sc["max_abs_ddg"] = sc.apply(lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)
    sc["tier_n"] = sc["comavi_tier"].astype(str).str.extract(r"(\d)").astype(float)

    L = sc[sc["label"].isin(["pathogenic", "benign"])].dropna(subset=["tier_n"]).copy()
    L["pos"] = L["label"].eq("pathogenic")
    L["tier_score"] = 5 - L["tier_n"]          # higher = stronger structural evidence
    return L


def main() -> int:
    L = load_labeled()
    z = lambda x: (x - x.mean()) / x.std(ddof=0)
    combined = z(L["tier_score"]) + z(L["max_abs_ddg"])

    a_tier = auc(L["tier_score"], L["pos"])
    a_ddg = auc(L["max_abs_ddg"], L["pos"])
    a_comb = auc(combined, L["pos"])

    rng = np.random.default_rng(SEED)
    systems = L["system"].values
    uniq = np.unique(systems)
    deltas, tiers, ddgs = [], [], []
    for _ in range(N_BOOT):
        draw = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.where(systems == s)[0] for s in draw])
        b = L.iloc[idx]
        if b["pos"].nunique() < 2:
            continue
        bc = z(b["tier_score"]) + z(b["max_abs_ddg"])
        deltas.append(auc(bc, b["pos"]) - auc(b["max_abs_ddg"], b["pos"]))
        tiers.append(auc(b["tier_score"], b["pos"]))
        ddgs.append(auc(b["max_abs_ddg"], b["pos"]))
    deltas, tiers, ddgs = map(np.asarray, (deltas, tiers, ddgs))

    def ci(a):
        return [round(float(np.percentile(a, 2.5)), 3), round(float(np.percentile(a, 97.5)), 3)]

    rho = stats.spearmanr(L["tier_n"], L["max_abs_ddg"])
    out = {
        "n": int(len(L)),
        "n_pathogenic": int(L["pos"].sum()),
        "n_benign": int((~L["pos"]).sum()),
        "n_systems": int(len(uniq)),
        "auc_tier_alone": round(a_tier, 3),
        "auc_energy_alone": round(a_ddg, 3),
        "auc_combined": round(a_comb, 3),
        "delta_auc_combined_minus_energy": round(a_comb - a_ddg, 4),
        "delta_auc_cluster_ci95": ci(deltas),
        "auc_tier_cluster_ci95": ci(tiers),
        "auc_energy_cluster_ci95": ci(ddgs),
        "spearman_tier_vs_energy": round(float(rho.statistic), 3),
        "spearman_p": float(f"{rho.pvalue:.3g}"),
        "n_bootstrap": int(len(deltas)),
        "note": (
            "Equal-weight z-sum, not a fitted combination: any weight fit on these "
            "49 variants and scored on the same 49 would be optimistically biased. "
            "Cluster bootstrap resamples whole systems."
        ),
    }
    dest = REPO / "reference_outputs/COMAVI_tier_auc_decomposition.json"
    dest.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    print(f"\nwrote {dest.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
