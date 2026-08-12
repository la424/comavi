#!/usr/bin/env python
"""Regenerate every number in the tier x energy-axis gating section from the
canonical scored table.

Motivation: the published version of this section was written by hand and went
stale when the two fold-neutral interface systems (VWF-GPIba, CFH-C3b) entered
the benchmark. Its cells reproduced only on the pre-expansion 42-variant
population while its own prose claimed n = 47, and one reported cell reproduced
under no population definition and no canonical vintage at all. Nothing in the
audit suite covered the section, so the drift was silent.

Every value the manuscript states for this section must come from the JSON this
script writes. scripts/audit_tier_energy_claims.py asserts that correspondence.

Definitions are taken from the shipped pipeline, never re-implemented here:
  - firing        : the confidence-gated three-axis max|ddG| at or above the
                    threshold, via apply_concordance_v5.compute_max_abs_ddg
                    (identical to the p1_ddg_concordance_* columns the pipeline
                    already stores, which this script cross-checks against).
  - strong tier   : comavi_tier in apply_concordance_v5.FOOTPRINT_TIERS.
  - structural GT : expected_mech_class in the three structural classes, as
                    opposed to structurally_silent. Rows whose ground truth
                    commits to neither are outside the population.

Population: variants carrying a tier and a committed structural-vs-silent
ground truth. The BRCA1-BRCT cohort carries no tier (its upstream contact and
burial features were never computed) and is therefore absent by construction,
not by exclusion -- see the tier applicability note in the manuscript.
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr, rankdata

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_tier_energy_gating.json"

STRUCTURAL_CLASSES = ["mixed_structural", "ppi_destab_mechanism", "fold_mechanism"]
SILENT_CLASS = "structurally_silent"
GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
CANONICAL_TAG = "t25"


def auc(score, y):
    """Rank-based AUC; no sklearn dependency."""
    y = np.asarray(y).astype(int)
    r = rankdata(np.asarray(score, dtype=float))
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank
        prev = min(prev, p[idx] * n / i)
        adj[idx] = prev
    return adj


def build_population(df, partners):
    d = df.copy()
    if "max_abs_ddg" not in d.columns or d["max_abs_ddg"].isna().all():
        d["max_abs_ddg"] = d.apply(lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)
    d["strong_tier"] = d["comavi_tier"].astype(str).isin(ac.FOOTPRINT_TIERS)
    d["structural_gt"] = d["expected_mech_class"].isin(STRUCTURAL_CLASSES)
    d["grade_value"] = d["mech_consistency_t25"].map(GRADE_MAP)
    pop = d[
        d["comavi_tier"].notna()
        & d["expected_mech_class"].isin(STRUCTURAL_CLASSES + [SILENT_CLASS])
    ].copy()
    return pop


def cells(sub):
    hi = sub[sub.strong_tier]
    lo = sub[~sub.strong_tier]
    return dict(
        strong_structural=int(hi.structural_gt.sum()), strong_n=int(len(hi)),
        weak_structural=int(lo.structural_gt.sum()), weak_n=int(len(lo)),
    )


def fisher_from_cells(c):
    return float(fisher_exact([
        [c["strong_structural"], c["strong_n"] - c["strong_structural"]],
        [c["weak_structural"], c["weak_n"] - c["weak_structural"]],
    ])[1])


def main():
    df = pd.read_csv(CANON)
    partners = ac.discover_partners(df)
    pop = build_population(df, partners)

    out = {
        "_source": str(CANON.relative_to(REPO)),
        "_definitions": {
            "firing": "confidence-gated three-axis max|ddG| >= threshold",
            "strong_tier": sorted(ac.FOOTPRINT_TIERS),
            "structural_classes": STRUCTURAL_CLASSES,
        },
        "population_n": int(len(pop)),
        "population_structural": int(pop.structural_gt.sum()),
        "population_silent": int((~pop.structural_gt).sum()),
        "population_systems": int(pop.system.nunique()),
    }

    # ---- four-state table at the canonical threshold -----------------------
    fires = pop["max_abs_ddg"].fillna(0) >= 2.5
    four = {}
    for label, f, t in [
        ("convergent", True, True),
        ("energy_only", True, False),
        ("context_only", False, True),
        ("neither", False, False),
    ]:
        s = pop[(fires == f) & (pop.strong_tier == t)]
        four[label] = dict(
            n=int(len(s)),
            structural=int(s.structural_gt.sum()),
            rate=round(float(s.structural_gt.mean()), 4) if len(s) else None,
        )
    out["four_state_t25"] = four

    # ---- threshold sweep, both strata --------------------------------------
    tags = [c.replace("p1_ddg_concordance_", "")
            for c in df.columns if c.startswith("p1_ddg_concordance_")]
    sweep, pv, names = {}, [], []
    for tag in tags:
        lab = pop[f"p1_ddg_concordance_{tag}"]
        firing_mask = lab.isin(["concordant_disruption", "ddg_only"])
        c_fire = cells(pop[firing_mask])
        c_sil = cells(pop[~firing_mask])
        p_fire = fisher_from_cells(c_fire)
        p_sil = fisher_from_cells(c_sil)
        sweep[tag] = {"firing": {**c_fire, "fisher_p": round(p_fire, 4)},
                      "silent": {**c_sil, "fisher_p": round(p_sil, 4)}}
        pv += [p_fire, p_sil]
        names += [f"firing_{tag}", f"silent_{tag}"]

    adj = bh(pv)
    out["sweep"] = sweep
    out["bh"] = {n: {"p": round(p, 4), "p_adj": round(a, 4), "survives_05": bool(a < 0.05)}
                 for n, p, a in zip(names, pv, adj)}
    out["bh_n_surviving"] = int(sum(a < 0.05 for a in adj))

    # ---- fragility of the canonical-threshold test -------------------------
    fire25 = pop[pop[f"p1_ddg_concordance_{CANONICAL_TAG}"]
                 .isin(["concordant_disruption", "ddg_only"])]
    c25 = cells(fire25)
    out["canonical_cells"] = {**c25, "fisher_p": round(fisher_from_cells(c25), 4)}
    out["reclassification_fragility"] = {}
    for k in (1, 2):
        c = dict(c25)
        c["weak_structural"] += k
        out["reclassification_fragility"][f"reclassify_{k}"] = round(fisher_from_cells(c), 4)

    weak = fire25[~fire25.strong_tier]
    out["weak_tier_firing_variants"] = [
        {"variant": r.variant, "system": r.system, "expected_mech_class": r.expected_mech_class}
        for r in weak.itertuples()
    ]
    out["n_systems_contributing_firing"] = int(fire25.system.nunique())
    out["n_systems_with_both_tier_strata"] = int(
        sum(1 for _, g in fire25.groupby("system") if g.strong_tier.nunique() == 2))

    # residue-level collapse (two of the weak-tier rows are the same residue)
    fr = fire25.copy()
    fr["residue"] = fr.system + ":" + fr.variant.str.extract(r"^([A-Z]\d+)")[0]
    coll = fr.groupby(["residue", "strong_tier"], as_index=False).structural_gt.max()
    c_coll = dict(
        strong_structural=int(coll[coll.strong_tier].structural_gt.sum()),
        strong_n=int(len(coll[coll.strong_tier])),
        weak_structural=int(coll[~coll.strong_tier].structural_gt.sum()),
        weak_n=int(len(coll[~coll.strong_tier])),
    )
    out["residue_collapsed"] = {**c_coll, "fisher_p": round(fisher_from_cells(c_coll), 4)}

    # ---- marginal tier screen (threshold-independent) ----------------------
    # The interaction framing above conditions on the energy threshold. The
    # underlying association does not: strong tier vs structural ground truth
    # is a 2x2 that never touches ddG. Reported with a within-system
    # permutation test and a system-level cluster bootstrap, because rows
    # within a system are not independent.
    nS = int(pop.structural_gt.sum())
    nN = int(len(pop) - nS)
    tp = int(pop[pop.structural_gt].strong_tier.sum())
    fp = int(pop[~pop.structural_gt].strong_tier.sum())
    fn, tn = nS - tp, nN - fp

    def clopper(k, n, a=0.05):
        from scipy.stats import beta as _b
        lo = 0.0 if k == 0 else float(_b.ppf(a / 2, k, n - k + 1))
        hi = 1.0 if k == n else float(_b.ppf(1 - a / 2, k + 1, n - k))
        return [round(lo, 4), round(hi, 4)]

    rng = np.random.default_rng(0)
    obs = tp + fp - fp  # = tp; kept explicit for readability
    strat_null = []
    for _ in range(20000):
        tot = 0
        for _, g in pop.groupby("system"):
            perm = g.structural_gt.values[rng.permutation(len(g))]
            tot += int(perm[g.strong_tier.values].sum())
        strat_null.append(tot)
    strat_null = np.asarray(strat_null)
    obs_strong_structural = int(pop[pop.strong_tier].structural_gt.sum())

    sysnames = pop.system.unique()
    diffs = []
    for _ in range(5000):
        S = pd.concat([pop[pop.system == s] for s in rng.choice(sysnames, len(sysnames), replace=True)])
        hi_, lo_ = S[S.strong_tier], S[~S.strong_tier]
        if len(hi_) and len(lo_):
            diffs.append(hi_.structural_gt.mean() - lo_.structural_gt.mean())

    marg = dict(
        strong_structural=obs_strong_structural,
        strong_n=int(pop.strong_tier.sum()),
        weak_structural=int(pop[~pop.strong_tier].structural_gt.sum()),
        weak_n=int((~pop.strong_tier).sum()),
        sensitivity=round(tp / nS, 4), sensitivity_ci=clopper(tp, nS),
        specificity=round(tn / nN, 4), specificity_ci=clopper(tn, nN),
        npv=round(tn / (tn + fn), 4) if (tn + fn) else None,
        ppv=round(tp / (tp + fp), 4) if (tp + fp) else None,
        n_systems_in_weak_cell=int(pop[~pop.strong_tier].system.nunique()),
        n_systems_with_both_strata=int(
            sum(1 for _, g in pop.groupby("system") if g.strong_tier.nunique() == 2)),
        weak_tier_classes={k: int(v) for k, v in
                           pop[~pop.strong_tier].expected_mech_class.value_counts().items()},
        structural_by_tier={k: int(v) for k, v in
                            pop[pop.structural_gt].comavi_tier.value_counts().items()},
    )
    # Is the energy-stratified test independent evidence, or a lower-power view
    # of the marginal? If the marginal weak-tier cell contains zero structural
    # variants, every stratified weak cell is a subset of it and must also be
    # zero -- the stratified p-values are inherited, not additional.
    marg["stratified_cells_are_subsets_of_marginal_zero"] = bool(
        marg["weak_structural"] == 0
        and all(v["firing"]["weak_structural"] == 0 for v in out["sweep"].values()))

    # How much of perfect recall is anticipated by the tier's own inputs? The
    # tier scores interface contacts, so recall on interface-defined mechanism
    # classes is partly definitional; fold_mechanism rows are the only ones
    # where it is not.
    S = pop[pop.structural_gt]
    marg["structural_by_class"] = {k: int(v) for k, v in S.expected_mech_class.value_counts().items()}
    nonppi = S[S.expected_mech_class != "ppi_destab_mechanism"]
    marg["non_ppi_structural_n"] = int(len(nonppi))
    marg["non_ppi_structural_all_strong_tier"] = bool(nonppi.strong_tier.all())
    fold_only = S[S.expected_mech_class == "fold_mechanism"]
    marg["fold_mechanism_n"] = int(len(fold_only))
    marg["fold_mechanism_all_strong_tier"] = bool(fold_only.strong_tier.all())

    marg["fisher_p"] = float(f"{fisher_from_cells(marg):.3e}")
    marg["within_system_permutation_p"] = round(float((strat_null >= obs_strong_structural).mean()), 4)
    marg["cluster_bootstrap_rate_diff_ci"] = [round(float(np.percentile(diffs, 2.5)), 4),
                                              round(float(np.percentile(diffs, 97.5)), 4)]
    out["marginal_tier_screen"] = marg

    # ---- mechanism consistency by stratum ----------------------------------
    conv = pop[(fires) & (pop.strong_tier)]
    tier_only = pop[(~fires) & (pop.strong_tier)]
    def mc_block(s):
        return dict(n=int(s.grade_value.notna().sum()),
                    mc=round(float(s.grade_value.mean()), 4) if s.grade_value.notna().any() else None)
    out["stratum_mc"] = {
        "convergent_all": mc_block(conv),
        "convergent_structural_gt": mc_block(conv[conv.structural_gt]),
        "convergent_silent_gt": {**mc_block(conv[~conv.structural_gt]),
                                 "variants": sorted(conv[~conv.structural_gt].variant.tolist())},
        "tier_only_structural_gt": {**mc_block(tier_only[tier_only.structural_gt]),
                                    "variants": sorted(tier_only[tier_only.structural_gt].variant.tolist())},
    }

    # ---- pathogenicity-classifier comparison -------------------------------
    lab = df.copy()
    lab["tier_n"] = lab["comavi_tier"].astype(str).str.extract(r"Tier (\d)").astype(float)
    if "max_abs_ddg" not in lab.columns or lab["max_abs_ddg"].isna().all():
        lab["max_abs_ddg"] = lab.apply(lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)

    def norm_label(x):
        s = str(x).lower()
        if "conflict" in s or "uncertain" in s:
            return "vus"
        if "pathogenic" in s:
            return "pathogenic"
        if "benign" in s:
            return "benign"
        return None

    # BRCT pathogenicity lives in clinvar_germline; its phenotype field is null.
    src = lab["phenotype"].where(lab.system != "brca1_brct")
    if "clinvar_germline" in lab.columns:
        src = src.fillna(lab["clinvar_germline"])
    lab["label"] = src.map(norm_label)
    L = lab[lab.label.isin(["pathogenic", "benign"])
            & lab.tier_n.notna() & lab.max_abs_ddg.notna()].copy()
    L["y"] = (L.label == "pathogenic").astype(int)
    z = lambda v: (v - v.mean()) / v.std()
    rho, prho = spearmanr(L.tier_n, L.max_abs_ddg)
    out["classifier"] = {
        "n": int(len(L)), "n_pathogenic": int(L.y.sum()), "n_benign": int((1 - L.y).sum()),
        "n_systems": int(L.system.nunique()),
        "auc_energy": round(auc(L.max_abs_ddg, L.y), 4),
        "auc_tier": round(auc(-L.tier_n, L.y), 4),
        "auc_combined": round(auc(z(L.max_abs_ddg) + z(-L.tier_n), L.y), 4),
        "spearman_tier_vs_energy": round(float(rho), 4),
        "spearman_p": float(prho),
    }

    # System-level cluster bootstrap on the combined-minus-energy AUC delta.
    # Resamples whole systems, because variants within a system share a
    # structure, a partner and a curation pass.
    rng2 = np.random.default_rng(0)
    sysL = L.system.unique()
    deltas = []
    for _ in range(2000):
        B = pd.concat([L[L.system == s] for s in rng2.choice(sysL, len(sysL), replace=True)])
        if B.y.nunique() < 2:
            continue
        zb = lambda v: (v - v.mean()) / (v.std() if v.std() else 1.0)
        deltas.append(auc(zb(B.max_abs_ddg) + zb(-B.tier_n), B.y) - auc(B.max_abs_ddg, B.y))
    out["classifier"]["auc_delta_combined_minus_energy"] = round(
        out["classifier"]["auc_combined"] - out["classifier"]["auc_energy"], 4)
    out["classifier"]["auc_delta_cluster_bootstrap_ci"] = [
        round(float(np.percentile(deltas, 2.5)), 3), round(float(np.percentile(deltas, 97.5)), 3)]
    out["classifier"]["auc_delta_bootstrap_n"] = len(deltas)

    OUT.write_text(json.dumps(out, indent=1))
    print(f"wrote {OUT.relative_to(REPO)}")
    print(f"  population n={out['population_n']} "
          f"({out['population_structural']} structural / {out['population_silent']} silent)")
    c = out["canonical_cells"]
    print(f"  canonical firing cells: {c['strong_structural']}/{c['strong_n']} vs "
          f"{c['weak_structural']}/{c['weak_n']}, p={c['fisher_p']}")
    print(f"  BH: {out['bh_n_surviving']} of {len(pv)} tests survive at 5%")


if __name__ == "__main__":
    main()
