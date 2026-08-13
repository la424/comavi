#!/usr/bin/env python
"""
verify_tier_construction.py — independent reproduction of the structural-context
tier from its published component formula, and the interface-bonus ablation.

WHY THIS EXISTS
The manuscript's circularity argument was originally written from mechanism-CLASS
labels: only 6 of the 17 structural-mechanism variants are binding-only class, so
it was asserted that the tier's interface bonus "accounts for 6 of the 17" and the
other 11 reach Tier 1-2 without it. That inference is invalid. Mechanism class is
ground truth about the variant; the interface bonus is awarded from pLDDT-gated
interface geometry, which many non-binding-class variants also have. The only
correct way to answer the question is to recompute the tier with the term removed.

Result: 15 of 17 structural variants receive a nonzero interface bonus, and
zeroing the term breaks perfect recall (17/17 -> 14/17) and the within-system
permutation test (p = 0.0095 -> 0.1294).

Also checks two claims that are ALGEBRAIC IDENTITIES rather than empirical tests,
and must not be reported as robustness evidence:
  (1) "requiring Tier 1-2 costs no sensitivity" — forced by the empty weak-tier
      structural cell, for ANY energy rule, at every threshold;
  (2) leave-one-system-out preserves the empty cell — removing rows can never
      create a weak-tier structural variant.

Regenerates reference_outputs/COMAVI_tier_construction.json.
"""
from __future__ import annotations

import json
import pathlib
import sys
from math import comb

import numpy as np
import pandas as pd
from scipy.stats import beta as _beta
from scipy.stats import fisher_exact

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs/scored_61var_canonical.csv"
LEDGER = REPO / "reference_outputs/COMAVI_evidence_ledger.csv"
OUT = REPO / "reference_outputs/COMAVI_tier_construction.json"

STRUCTURAL_CLASSES = ["mixed_structural", "ppi_destab_mechanism", "fold_mechanism"]
SILENT_CLASS = "structurally_silent"
FIRING_LABELS = ["concordant_disruption", "ddg_only"]

# component constants, mirrored from scripts/comavi_v7/config.py + mechanism.py
TIER_THRESHOLDS = [(5.0, "Tier 1"), (3.0, "Tier 2"), (1.5, "Tier 3")]
DEFAULT_TIER = "Tier 4"
DISRUPTION_POINTS = [(20, 4.0), (10, 3.0), (4, 2.0), (1, 1.0)]
BURIAL_RANK = {"buried_core": 2, "partially_buried": 1, "surface_exposed": 0}


def _f(v, d=0.0):
    return d if v is None or pd.isna(v) else float(v)


def tier_score(row, partners, use_interface_bonus=True):
    """Reimplementation of comavi_v7.mechanism.compute_disruption_score."""
    score = 0.0
    sev = _f(row.get("substitution_severity"))
    mono_c = _f(row.get("monomer_n_contacts"))

    max_inter, gated_iface = 0.0, []
    for pl in partners:
        if _f(row.get(f"multi_{pl}_plddt")) < 50:
            continue
        iv = row.get(f"multi_{pl}_is_interface")
        if iv is not None and not pd.isna(iv) and bool(iv) and str(iv).lower() not in ("0", "0.0", "false"):
            gated_iface.append(pl)
        ic = row.get(f"multi_{pl}_inter_contacts")
        if pd.notna(ic):
            max_inter = max(max_inter, float(ic))

    disruption = round(sev * (mono_c + max_inter), 2)
    for thresh, pts in DISRUPTION_POINTS:
        if disruption >= thresh:
            score += pts
            break

    if use_interface_bonus:
        n = len(gated_iface)
        score += 2.0 if n >= 2 else (1.5 if n == 1 else 0.0)

    best_rank = BURIAL_RANK.get(str(row.get("monomer_burial")), 0)
    for pl in partners:
        b = row.get(f"multi_{pl}_burial")
        if pd.notna(b) and _f(row.get(f"multi_{pl}_plddt")) >= 50:
            best_rank = max(best_rank, BURIAL_RANK.get(str(b), 0))
    score += 2.0 if best_rank == 2 else (1.0 if best_rank == 1 else 0.0)

    bp = row.get("best_plddt")
    if pd.notna(bp):
        if float(bp) < 50:
            score *= 0.4
        elif float(bp) < 70:
            score *= 0.7

    return round(score, 2), len(gated_iface)


def to_tier(s):
    for t, lab in TIER_THRESHOLDS:
        if s >= t:
            return lab
    return DEFAULT_TIER


def cp_interval(k, n, alpha=0.05):
    lo = 0.0 if k == 0 else float(_beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(_beta.ppf(1 - alpha / 2, k + 1, n - k))
    return [round(lo, 3), round(hi, 3)]


def exact_within_system_perm(pop, strong_col):
    """Exact (enumerated, not Monte Carlo) within-system permutation p-value.

    Under H0 the structural labels are exchangeable WITHIN each system, so the
    count landing in the strong-tier cell is a sum of independent hypergeometrics.
    Convolve the per-system distributions and take the upper tail.
    """
    obs = int((pop[strong_col] & pop["structural_gt"]).sum())
    dists = []
    for _, g in pop.groupby("system"):
        n, k, ns = len(g), int(g["structural_gt"].sum()), int(g[strong_col].sum())
        d = {}
        for x in range(max(0, k - (n - ns)), min(k, ns) + 1):
            d[x] = comb(ns, x) * comb(n - ns, k - x) / comb(n, k)
        dists.append(d)
    tot = {0: 1.0}
    for d in dists:
        nt = {}
        for s0, p0 in tot.items():
            for x, px in d.items():
                nt[s0 + x] = nt.get(s0 + x, 0.0) + p0 * px
        tot = nt
    return float(sum(p for s, p in tot.items() if s >= obs)), obs


def screen_stats(pop, strong_col):
    S = pop["structural_gt"]
    a = int((pop[strong_col] & S).sum())
    b = int((~pop[strong_col] & S).sum())
    c = int((pop[strong_col] & ~S).sum())
    d = int((~pop[strong_col] & ~S).sum())
    p, _ = exact_within_system_perm(pop, strong_col)
    return {
        "strong_structural": a, "weak_structural": b,
        "strong_silent": c, "weak_silent": d,
        "sensitivity": round(a / (a + b), 3), "sensitivity_ci": cp_interval(a, a + b),
        "specificity": round(d / (c + d), 3), "specificity_ci": cp_interval(d, c + d),
        # 6 dp, not 3 sig-figs: the audit formats this to 4 dp, and storing
        # 3 sig-figs first would double-round (0.006546 -> 0.00655 -> 0.0066).
        "fisher_p": round(float(fisher_exact([[a, c], [b, d]])[1]), 6),
        "within_system_permutation_p": round(p, 4),
    }


def auc(score, y):
    y = np.asarray(y, dtype=bool)
    s = np.asarray(score, dtype=float)
    pos, neg = s[y], s[~y]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def main():
    df = pd.read_csv(CANON)
    partners = ac.discover_partners(df)

    d = df.copy()
    d["max_abs_ddg"] = d.apply(lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)
    d["strong_tier"] = d["comavi_tier"].astype(str).isin(ac.FOOTPRINT_TIERS)
    d["structural_gt"] = d["expected_mech_class"].isin(STRUCTURAL_CLASSES)
    pop = d[d["comavi_tier"].notna()
            & d["expected_mech_class"].isin(STRUCTURAL_CLASSES + [SILENT_CLASS])].copy()

    rec = pop.apply(lambda r: tier_score(r, partners, True), axis=1)
    pop["tier_recomputed"] = [to_tier(s) for s, _ in rec]
    pop["n_iface_bonus"] = [n for _, n in rec]
    fidelity = int((pop["tier_recomputed"] == pop["comavi_tier"]).sum())
    assert fidelity == len(pop), (
        f"tier reimplementation does not reproduce shipped comavi_tier "
        f"({fidelity}/{len(pop)}) — component formula has drifted from mechanism.py"
    )

    rec0 = pop.apply(lambda r: tier_score(r, partners, False), axis=1)
    pop["tier_no_iface"] = [to_tier(s) for s, _ in rec0]
    pop["strong_no_iface"] = pop["tier_no_iface"].isin(list(ac.FOOTPRINT_TIERS))

    S = pop[pop["structural_gt"]]
    out = {
        "population_n": len(pop),
        "structural_n": int(pop["structural_gt"].sum()),
        "silent_n": int((~pop["structural_gt"]).sum()),
        "reimplementation_fidelity": f"{fidelity}/{len(pop)}",
        "interface_bonus_reach": {
            "structural_with_bonus": int((S["n_iface_bonus"] >= 1).sum()),
            "structural_total": len(S),
            "by_class": {k: f"{int((g['n_iface_bonus'] >= 1).sum())}/{len(g)}"
                         for k, g in S.groupby("expected_mech_class")},
            "structural_without_bonus": S.loc[S["n_iface_bonus"] == 0, "variant"].tolist(),
            "note": ("15/17 — NOT the 6 binding-only-class variants. The manuscript's "
                     "original 6-of-17 claim inferred bonus reach from mechanism class "
                     "and was wrong."),
        },
        "screen_full_tier": screen_stats(pop, "strong_tier"),
        "screen_interface_bonus_zeroed": screen_stats(pop, "strong_no_iface"),
        "variants_demoted_by_ablation": pop.loc[
            pop["strong_tier"] & ~pop["strong_no_iface"],
            ["variant", "system", "expected_mech_class", "comavi_tier", "tier_no_iface"],
        ].to_dict("records"),
    }

    # ---- identity check 1: sensitivity invariance is forced, not observed -----
    forced = {}
    for tag, _spec in ac.THRESHOLD_SPECS:
        fires = pop[f"p1_ddg_concordance_{tag}"].isin(FIRING_LABELS)
        forced[tag] = int((fires & ~pop["strong_tier"] & pop["structural_gt"]).sum())
    out["sensitivity_invariance_is_identity"] = {
        "structural_variants_removed_by_requiring_strong_tier": forced,
        "note": ("Zero at every threshold BECAUSE the marginal weak-tier structural "
                 "cell is empty. Conjoining strong_tier to ANY energy rule removes "
                 "zero true positives. Do not report as 'no sensitivity cost'."),
    }

    # ---- identity check 2: LOSO cannot break the empty cell ------------------
    loso = {}
    for s in sorted(pop["system"].unique()):
        sub = pop[pop["system"] != s]
        loso[s] = {
            "n": len(sub),
            "weak_tier_structural": int((~sub["strong_tier"] & sub["structural_gt"]).sum()),
            "sensitivity": round(float(sub.loc[sub["structural_gt"], "strong_tier"].mean()), 3),
        }
    out["loso_is_identity"] = {
        "per_system": loso,
        "note": ("Row removal cannot CREATE a weak-tier structural variant when the "
                 "full-cohort cell is already 0. Valid as a check on the tier "
                 "GRADIENT (refit), not on the screen's perfect recall."),
    }

    # ---- specificity gain: where it lives -----------------------------------
    rng = np.random.default_rng(0)
    systems = pop["system"].unique()
    gain_rows = {}
    for tag, _spec in ac.THRESHOLD_SPECS:
        fires = pop[f"p1_ddg_concordance_{tag}"].isin(FIRING_LABELS)
        eo = pop[fires & ~pop["strong_tier"]]
        sil = pop[~pop["structural_gt"]]
        spec_e = float((~fires[sil.index]).mean())
        spec_c = float((~(fires[sil.index] & sil["strong_tier"])).mean())
        gain_rows[tag] = {
            "specificity_energy_only": round(spec_e, 3),
            "specificity_convergent": round(spec_c, 3),
            "gain": round(spec_c - spec_e, 3),
            "energy_only_cell_n": len(eo),
            "energy_only_systems": int(eo["system"].nunique()),
            "energy_only_distinct_residues": int(eo["variant"].str[1:-1].nunique()),
            "energy_only_variants": eo["variant"].tolist(),
        }
    ref = pop.assign(fires=pop["p1_ddg_concordance_t25"].isin(FIRING_LABELS))
    gains = []
    for _ in range(4000):
        pick = rng.choice(systems, len(systems), replace=True)
        R = pd.concat([ref[ref["system"] == s] for s in pick], ignore_index=True)
        sil = R[~R["structural_gt"]]
        if len(sil) == 0:
            continue
        spec_energy = float((~sil["fires"]).mean())
        spec_convergent = float((~(sil["fires"] & sil["strong_tier"])).mean())
        gains.append(spec_convergent - spec_energy)
    gains = np.array(gains)
    out["specificity_gain"] = {
        "by_threshold": gain_rows,
        "reference_t25_cluster_bootstrap_ci": [round(float(np.percentile(gains, 2.5)), 3),
                                              round(float(np.percentile(gains, 97.5)), 3)],
        "fraction_resamples_no_gain": round(float((gains <= 0).mean()), 3),
    }

    # ---- evidence-restricted views -----------------------------------------
    led = pd.read_csv(LEDGER)
    strict = {"E1_quantitative_energetic", "E2_quantitative_functional", "E3_qualitative_experimental"}
    best = (led.groupby(["system", "variant"])
            .agg(has_e13=("evidence_type", lambda s: bool(set(s) & strict)),
                 direct=("evidence_directness", lambda s: bool({"direct", "coupled"} & set(s))))
            .reset_index())
    pv = pop.merge(best, on=["system", "variant"], how="left")
    ev = {}
    for lab, mask in [("all", pd.Series(True, index=pv.index)),
                      ("E1_E3_only", pv["has_e13"].fillna(False)),
                      ("direct_or_coupled_only", pv["direct"].fillna(False))]:
        sub = pv[mask]
        a = int((sub["strong_tier"] & sub["structural_gt"]).sum())
        b = int((~sub["strong_tier"] & sub["structural_gt"]).sum())
        c = int((sub["strong_tier"] & ~sub["structural_gt"]).sum())
        dd = int((~sub["strong_tier"] & ~sub["structural_gt"]).sum())
        ev[lab] = {"n": len(sub), "sensitivity": round(a / (a + b), 3),
                   "specificity": round(dd / (c + dd), 3), "weak_tier_structural": b,
                   "fisher_p": float(f"{fisher_exact([[a, c], [b, dd]])[1]:.3g}")}
    out["evidence_restricted_views"] = ev

    # ---- same-cohort structural-mechanism detection AUCs --------------------
    tiern = {"Tier 1": 1, "Tier 2": 2, "Tier 3": 3, "Tier 4": 4}
    e = pop["max_abs_ddg"].fillna(0).to_numpy(float)
    t = (5 - pop["comavi_tier"].map(tiern)).to_numpy(float)
    y = pop["structural_gt"].to_numpy(bool)
    z = lambda v: (v - v.mean()) / v.std(ddof=0)  # noqa: E731
    sysarr = pop["system"].to_numpy()
    de, dt = [], []
    for _ in range(4000):
        pick = rng.choice(np.unique(sysarr), len(np.unique(sysarr)), replace=True)
        idx = np.concatenate([np.where(sysarr == s)[0] for s in pick])
        ys = y[idx]
        if ys.all() or not ys.any():
            continue
        ee, tt = e[idx], t[idx]
        cc = z(ee) + z(tt)
        de.append(auc(cc, ys) - auc(ee, ys))
        dt.append(auc(cc, ys) - auc(tt, ys))
    de, dt = np.array(de), np.array(dt)
    out["detection_auc"] = {
        "energy": round(auc(e, y), 3),
        "tier": round(auc(t, y), 3),
        "combined_zsum": round(auc(z(e) + z(t), y), 3),
        "delta_combined_minus_energy_ci": [round(float(np.percentile(de, 2.5)), 3),
                                           round(float(np.percentile(de, 97.5)), 3)],
        "delta_combined_minus_tier_ci": [round(float(np.percentile(dt, 2.5)), 3),
                                         round(float(np.percentile(dt, 97.5)), 3)],
        "note": ("Same-cohort, not externally validated. The tier's interface bonus "
                 "is partly definitional for structural class (15/17), so the "
                 "combination gain is not independent of ground truth."),
    }

    OUT.write_text(json.dumps(out, indent=1) + "\n")

    print(f"tier reimplementation fidelity: {out['reimplementation_fidelity']}")
    ib = out["interface_bonus_reach"]
    print(f"interface bonus reaches {ib['structural_with_bonus']}/{ib['structural_total']} "
          f"structural variants  {ib['by_class']}")
    for k in ("screen_full_tier", "screen_interface_bonus_zeroed"):
        s = out[k]
        print(f"{k:32} sens={s['sensitivity']} ({s['strong_structural']}/"
              f"{s['strong_structural'] + s['weak_structural']}) spec={s['specificity']} "
              f"fisher={s['fisher_p']} perm={s['within_system_permutation_p']}")
    print("demoted by ablation:",
          [f"{r['variant']}({r['system']})" for r in out["variants_demoted_by_ablation"]])
    print("sensitivity removed by requiring strong tier:",
          out["sensitivity_invariance_is_identity"]["structural_variants_removed_by_requiring_strong_tier"])
    g = out["specificity_gain"]["by_threshold"]["t25"]
    print(f"t25 specificity gain {g['gain']:+} from an energy-only cell of "
          f"{g['energy_only_cell_n']} variants / {g['energy_only_systems']} systems / "
          f"{g['energy_only_distinct_residues']} residues; cluster CI "
          f"{out['specificity_gain']['reference_t25_cluster_bootstrap_ci']}")
    a = out["detection_auc"]
    print(f"detection AUC energy={a['energy']} tier={a['tier']} combined={a['combined_zsum']} "
          f"delta-vs-energy CI {a['delta_combined_minus_energy_ci']}")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
