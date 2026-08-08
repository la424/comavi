#!/usr/bin/env python
"""Regenerate every number in the audited threshold / tier sections of the COMAVI paper.

This script exists because three claims drafted earlier did not survive audit and
were replaced. It regenerates the replacements, and it also regenerates the
falsification checks themselves so a reader can confirm the retractions.

RETRACTED AND REPLACED
  R1. "A benchmark weighted the other way would place the optimum lower."
      -> True but vacuous as stated. Replaced by the explicit crossover weight.
  R2. "The threshold-labile variants split perfectly by ground-truth class."
      -> Tautological: raising a threshold can only remove firings, so the
         direction each variant moves is DETERMINED by its class. Replaced by
         the non-forced part: which variants are labile, and their effect sizes.
  R3. "The canonical threshold corresponds to 3-4 kcal/mol of real destabilization."
      -> Not defensible. Both shipped fits regress FoldX on measured; reading a
         threshold back is a reverse prediction whose Fieller intervals are far
         too wide (instability index g > 0.1 for both fits). Replaced by a
         model-free count of what each threshold catches and misses.

Outputs land in reference_outputs/ as COMAVI_audit_*.{csv,json}.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, linregress, mannwhitneyu, spearmanr, t as tdist

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

GMAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
STRUCT = {"ppi_destab_mechanism", "mixed_structural", "fold_mechanism"}
REAL_DESTAB_KCAL = 1.0  # model-free definition of a real destabilizer
N_BOOT = 10000
SEED = 42


def load():
    df = pd.read_csv(REPO / "reference_outputs/scored_61var_canonical.csv")
    partners = [p for p in ac.discover_partners(df) if f"ddg_{p}_confident" in df.columns]
    tags = [t for t, _ in ac.THRESHOLD_SPECS]
    unobs = set(ac.unobservable_variants())
    g = df[~df.variant.isin(unobs)].copy()
    for t in tags:
        g[f"g_{t}"] = g[f"mech_consistency_{t}"].map(GMAP)
    com = g[g.expected_mech_class.isin(STRUCT | {"structurally_silent"})].copy()
    com["truth_struct"] = com.expected_mech_class.isin(STRUCT)
    com["max_abs_ddg"] = com.apply(lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)
    return df, g, com, partners, tags


def ddg_axis_fire(row, partners, spec):
    """Per-axis destabilizing fire under the shipped confidence gates.

    Matches scripts/build_tier_concordance.py: signed (destabilizing-only) ddG,
    not absolute value.
    """
    if isinstance(spec, dict):
        tm, tf, tb = spec["monomer"], spec["fold"], spec["binding"]
    else:
        tm = tf = tb = float(spec)
    mono = fold = bind = False
    mv = row.get("ddg_monomer")
    if pd.notna(mv) and bool(row.get("ddg_monomer_confident", False)):
        mono = float(mv) >= tm
    for p in partners:
        pl = str(p).lower()
        if not bool(row.get(f"ddg_{pl}_confident", False)):
            continue
        fv = row.get(f"ddg_fold_{pl}")
        if pd.notna(fv) and float(fv) >= tf:
            fold = True
        if bool(row.get(f"ddg_binding_{pl}_indistinguishable", False)):
            continue
        bv = row.get(f"ddg_binding_{pl}")
        if pd.notna(bv) and float(bv) >= tb:
            bind = True
    return mono or fold or bind


def operating_points(com, tags):
    gc = [f"g_{t}" for t in tags]
    recall = np.array([com[com.truth_struct][c].mean() for c in gc])
    reject = np.array([com[~com.truth_struct][c].mean() for c in gc])
    pooled = np.array([com[c].mean() for c in gc])
    return recall, reject, pooled, gc


def r1_crossover(recall, reject, tags):
    """R1 replacement: the exact structural-class weight at which the optimum moves."""
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if tags[int((mid * recall + (1 - mid) * reject).argmax())] == "t25":
            lo = mid
        else:
            hi = mid
    return hi


def r2_lability(com, gc):
    """R2 replacement: lability is forced in DIRECTION but not in membership or magnitude."""
    com = com.copy()
    com["n_distinct_grades"] = com[gc].nunique(axis=1)
    com["labile"] = com.n_distinct_grades > 1
    mono = int(sum(
        all(a >= b for a, b in zip(v, v[1:])) or all(a <= b for a, b in zip(v, v[1:]))
        for v in com[gc].values
    ))
    ls, lz = com[com.labile & com.truth_struct], com[com.labile & ~com.truth_struct]
    u, p = mannwhitneyu(ls.max_abs_ddg, lz.max_abs_ddg)
    return com, {
        "n_committed": int(len(com)),
        "n_stable": int((~com.labile).sum()),
        "n_labile": int(com.labile.sum()),
        "grades_monotone_in_threshold": mono,
        "direction_is_forced": bool(mono == len(com)),
        "labile_structural_n": int(len(ls)),
        "labile_silent_n": int(len(lz)),
        "labile_structural_median_ddg": round(float(ls.max_abs_ddg.median()), 3),
        "labile_silent_median_ddg": round(float(lz.max_abs_ddg.median()), 3),
        "stable_structural_median_ddg": round(float(com[~com.labile & com.truth_struct].max_abs_ddg.median()), 3),
        "stable_silent_median_ddg": round(float(com[~com.labile & ~com.truth_struct].max_abs_ddg.median()), 3),
        "mannwhitney_U": float(u),
        "mannwhitney_p": round(float(p), 5),
    }


def _fieller(x, y, target, alpha=0.05):
    n = len(x)
    f = linregress(x, y)
    xb = x.mean()
    sxx = ((x - xb) ** 2).sum()
    resid = y - (f.intercept + f.slope * x)
    s2 = (resid ** 2).sum() / (n - 2)
    tc = tdist.ppf(1 - alpha / 2, n - 2)
    point = (target - f.intercept) / f.slope
    gidx = (tc ** 2 * s2) / (f.slope ** 2 * sxx)
    if gidx >= 1:
        return point, np.nan, np.nan, gidx, f
    c = point - xb
    half = (tc * np.sqrt(s2) / f.slope) * np.sqrt(c ** 2 / sxx + (1 - gidx) / n) / (1 - gidx)
    ctr = xb + c / (1 - gidx)
    return point, ctr - half, ctr + half, gidx, f


def r3_calibration_and_modelfree():
    """R3: show the inversion is unstable, then replace it with model-free counting."""
    cal = pd.read_csv(REPO / "reference_outputs/COMAVI_delta_calibration_points.csv")
    d = cal.dropna(subset=["measured_kcal", "foldx_ddg"]).copy()
    d["abs_measured"] = d.measured_kcal.abs()
    d["abs_foldx"] = d.foldx_ddg.abs()
    fits = {
        "benchmark": d[d.in_benchmark == True],          # noqa: E712
        "pooled": cal[cal.g73 != True].dropna(subset=["measured_kcal", "foldx_ddg"]),  # noqa: E712
    }
    inv_rows = []
    for name, sub in fits.items():
        fwd = linregress(sub.measured_kcal, sub.foldx_ddg)
        rev = linregress(sub.foldx_ddg, sub.measured_kcal)
        for tag, spec in ac.THRESHOLD_SPECS:
            tv = spec["fold"] if isinstance(spec, dict) else spec
            pt, lo, hi, gidx, _ = _fieller(sub.measured_kcal.values, sub.foldx_ddg.values, tv)
            inv_rows.append(dict(
                fit=name, n=len(sub), r2=round(fwd.rvalue ** 2, 3), threshold=tag, foldx_t=tv,
                inverse_estimate=round(pt, 2), fieller_lo=round(lo, 2), fieller_hi=round(hi, 2),
                instability_g=round(gidx, 3),
                direct_regression_estimate=round(rev.slope * tv + rev.intercept, 2),
            ))
    inv = pd.DataFrame(inv_rows)

    real = d[d.abs_measured >= REAL_DESTAB_KCAL]
    mf_rows = []
    for tag, spec in ac.THRESHOLD_SPECS:
        tv = spec["fold"] if isinstance(spec, dict) else spec
        caught = real[real.abs_foldx >= tv]
        missed = real[real.abs_foldx < tv]
        biggest = missed.nlargest(1, "abs_measured")
        mf_rows.append(dict(
            threshold=tag, foldx_t=tv, n_real_destabilizers=len(real),
            caught=len(caught), missed=len(missed),
            measured_recall=round(len(caught) / len(real), 3),
            median_measured_caught=round(float(caught.abs_measured.median()), 2) if len(caught) else None,
            median_measured_missed=round(float(missed.abs_measured.median()), 2) if len(missed) else None,
            max_measured_missed=round(float(missed.abs_measured.max()), 2) if len(missed) else None,
            largest_miss=str(biggest.variant.iloc[0]) if len(biggest) else None,
            largest_miss_measured=round(float(biggest.measured_kcal.iloc[0]), 2) if len(biggest) else None,
            largest_miss_foldx=round(float(biggest.foldx_ddg.iloc[0]), 2) if len(biggest) else None,
        ))
    return inv, pd.DataFrame(mf_rows), len(real), len(d)


def gate_sweep(df, partners):
    base = df[df.comavi_tier.astype(str).str.startswith("Tier")
              & df.expected_mech_class.isin(STRUCT | {"structurally_silent"})].copy()
    base["tier_strong"] = base.comavi_tier.str.split().str[-1].astype(int) <= 2
    base["structural"] = base.expected_mech_class.isin(STRUCT)
    base["residue"] = base.variant.str.extract(r"^([A-Z]\d+)")[0]
    rows = []
    for tag, spec in ac.THRESHOLD_SPECS:
        b = base.copy()
        b["fires"] = b.apply(lambda r: ddg_axis_fire(r, partners, spec), axis=1)
        rec = {"threshold": tag}
        for key, sub in [("firing", b[b.fires]), ("silent", b[~b.fires])]:
            A = int((sub.tier_strong & sub.structural).sum())
            B = int((sub.tier_strong & ~sub.structural).sum())
            C = int((~sub.tier_strong & sub.structural).sum())
            D = int((~sub.tier_strong & ~sub.structural).sum())
            rec.update({
                f"{key}_strong_structural": A, f"{key}_strong_n": A + B,
                f"{key}_weak_structural": C, f"{key}_weak_n": C + D,
                f"{key}_fisher_p": round(float(fisher_exact([[A, B], [C, D]])[1]), 4),
            })
        wf = b[b.fires & ~b.tier_strong]
        rec["weak_firing_n_variants"] = int(len(wf))
        rec["weak_firing_n_residues"] = int(wf.residue.nunique())
        rec["weak_firing_n_systems"] = int(wf.system.nunique())
        rows.append(rec)
    sweep = pd.DataFrame(rows)

    b25 = base.copy()
    b25["fires"] = b25.apply(lambda r: ddg_axis_fire(r, partners, 2.5), axis=1)
    f25 = b25[b25.fires]
    coll = f25.groupby(["system", "residue"]).agg(
        tier_strong=("tier_strong", "first"), structural=("structural", "max")).reset_index()
    A = int((coll.tier_strong & coll.structural).sum())
    B = int((coll.tier_strong & ~coll.structural).sum())
    C = int((~coll.tier_strong & coll.structural).sum())
    D = int((~coll.tier_strong & ~coll.structural).sum())
    frag = {
        "weak_cell_n_variants": int((~f25.tier_strong).sum()),
        "weak_cell_n_residues": int(f25[~f25.tier_strong].residue.nunique()),
        "weak_cell_n_systems": int(f25[~f25.tier_strong].system.nunique()),
        "weak_cell_variants": sorted(f25[~f25.tier_strong].variant.tolist()),
        "residue_collapsed_2x2": [[A, B], [C, D]],
        "residue_collapsed_fisher_p": round(float(fisher_exact([[A, B], [C, D]])[1]), 4),
        "systems_with_both_tier_strata_among_firing": int(
            (f25.groupby("system").tier_strong.nunique() > 1).sum()),
        "n_systems_among_firing": int(f25.system.nunique()),
    }
    for flip in (1, 2, 3):
        A0 = A if False else int((f25.tier_strong & f25.structural).sum())
        B0 = int((f25.tier_strong & ~f25.structural).sum())
        C0 = int((~f25.tier_strong & f25.structural).sum())
        D0 = int((~f25.tier_strong & ~f25.structural).sum())
        if flip <= C0 + D0:
            frag[f"fisher_p_if_{flip}_weak_were_structural"] = round(
                float(fisher_exact([[A0, B0], [C0 + flip, D0 - flip]])[1]), 4)
    return sweep, frag


def tier_gradient(df):
    ta = df[df.comavi_tier.astype(str).str.startswith("Tier")].copy()
    ta["tier_n"] = ta.comavi_tier.str.split().str[-1].astype(int)

    def norm(v):
        s = str(v).lower()
        if "conflict" in s or "uncertain" in s or "vus" in s:
            return "vus"
        if "pathogenic" in s:
            return "pathogenic"
        if "benign" in s:
            return "benign"
        return None

    ta["plab"] = ta.phenotype.map(norm)
    pb = ta[ta.plab.isin(["pathogenic", "benign"])].copy()
    pb["is_path"] = (pb.plab == "pathogenic").astype(int)
    rho, p = spearmanr(pb.tier_n, pb.is_path)
    loso = []
    for s in sorted(pb.system.unique()):
        d2 = pb[pb.system != s]
        if d2.tier_n.nunique() < 2:
            continue
        r2, p2 = spearmanr(d2.tier_n, d2.is_path)
        loso.append(dict(dropped_system=s, n=len(d2), rho=round(float(r2), 3), p=round(float(p2), 4)))
    loso = pd.DataFrame(loso)
    rng = np.random.default_rng(11)
    us = pb.system.unique()
    boot = []
    for _ in range(5000):
        pick = rng.choice(us, len(us), replace=True)
        idx = np.concatenate([np.where(pb.system.values == s)[0] for s in pick])
        d2 = pb.iloc[idx]
        if d2.tier_n.nunique() < 2:
            continue
        boot.append(spearmanr(d2.tier_n, d2.is_path)[0])
    boot = np.array(boot)
    grad = [dict(tier=int(t), pathogenic=int(pb[pb.tier_n == t].is_path.sum()),
                 n=int((pb.tier_n == t).sum()),
                 frac=round(float(pb[pb.tier_n == t].is_path.mean()), 3))
            for t in sorted(pb.tier_n.unique())]
    return {
        "n": int(len(pb)), "spearman_rho": round(float(rho), 4), "spearman_p": round(float(p), 5),
        "gradient": grad,
        "loso_worst_p": float(loso.p.max()), "loso_worst_system": str(loso.loc[loso.p.idxmax(), "dropped_system"]),
        "loso_all_significant": bool((loso.p < 0.05).all()),
        "cluster_boot_rho_mean": round(float(np.nanmean(boot)), 3),
        "cluster_boot_ci": [round(float(np.nanpercentile(boot, 2.5)), 3),
                            round(float(np.nanpercentile(boot, 97.5)), 3)],
        "cluster_boot_p_rho_negative": round(float(np.mean(boot < 0)), 4),
    }, loso


def bootstrap_tradeoff(com, gc):
    rng = np.random.default_rng(SEED)
    S = com[com.truth_struct][gc].values
    Z = com[~com.truth_struct][gc].values
    bs = np.array([S[rng.integers(0, len(S), len(S))].mean(0) for _ in range(N_BOOT)])
    bz = np.array([Z[rng.integers(0, len(Z), len(Z))].mean(0) for _ in range(N_BOOT)])
    sysid = com.system.values
    uniq = np.unique(sysid)
    V = com[gc].values

    def clus(mask, B=4000):
        out = []
        for _ in range(B):
            pick = rng.choice(uniq, len(uniq), replace=True)
            parts = [np.where((sysid == s) & mask)[0] for s in pick]
            idx = np.concatenate([p for p in parts if len(p)])
            if len(idx):
                out.append(V[idx].mean(0))
        return np.array(out)

    cs, cz = clus(com.truth_struct.values), clus(~com.truth_struct.values)
    return {
        "variant_boot_recall_ci": [[round(float(np.percentile(bs[:, i], 2.5)), 3),
                                    round(float(np.percentile(bs[:, i], 97.5)), 3)] for i in range(bs.shape[1])],
        "variant_boot_reject_ci": [[round(float(np.percentile(bz[:, i], 2.5)), 3),
                                    round(float(np.percentile(bz[:, i], 97.5)), 3)] for i in range(bz.shape[1])],
        "p_recall_drops": round(float((bs[:, 0] > bs[:, -1]).mean()), 4),
        "p_reject_rises": round(float((bz[:, -1] > bz[:, 0]).mean()), 4),
        "p_curves_cross": round(float(((bs[:, 0] > bz[:, 0]) & (bz[:, -1] > bs[:, -1])).mean()), 4),
        "cluster_p_recall_drops": round(float((cs[:, 0] > cs[:, -1]).mean()), 4),
        "cluster_p_reject_rises": round(float((cz[:, -1] > cz[:, 0]).mean()), 4),
        "cluster_recall_ci_t10": [round(float(np.percentile(cs[:, 0], 2.5)), 3),
                                  round(float(np.percentile(cs[:, 0], 97.5)), 3)],
        "cluster_reject_ci_tSAP": [round(float(np.percentile(cz[:, -1], 2.5)), 3),
                                   round(float(np.percentile(cz[:, -1], 97.5)), 3)],
        "note": "Per-variant grades are monotone in threshold by construction, so the "
                "bootstrap confirms the MAGNITUDE of the tradeoff, not its shape.",
    }


def main():
    out = REPO / "reference_outputs"
    df, g, com, partners, tags = load()
    recall, reject, pooled, gc = operating_points(com, tags)

    ops = pd.DataFrame({
        "threshold": tags,
        "foldx_t": [s["fold"] if isinstance(s, dict) else s for _, s in ac.THRESHOLD_SPECS],
        "recall_structural": recall.round(4),
        "rejection_silent": reject.round(4),
        "pooled_MC": pooled.round(4),
        "balanced_MC": ((recall + reject) / 2).round(4),
    })
    com2, r2 = r2_lability(com, gc)
    inv, mf, n_real, n_cal = r3_calibration_and_modelfree()
    sweep, frag = gate_sweep(df, partners)
    grad, loso = tier_gradient(df)
    boot = bootstrap_tradeoff(com, gc)

    summary = {
        "denominators": {
            "n_rows": int(len(df)), "n_unobservable": int(len(df) - len(g)),
            "n_graded": int(g[f"g_{tags[3]}"].notna().sum()),
            "n_committed": int(len(com)),
            "n_structural": int(com.truth_struct.sum()),
            "n_silent": int((~com.truth_struct).sum()),
            "note": "graded == committed == 57; the 2 class-uncommitted variants carry no grade, "
                    "so every headline denominator coincides.",
        },
        "headline_MC_t25": round(float(com[f"g_{tags[3]}"].mean()), 4),
        "threshold_averaged_MC": round(float(pooled.mean()), 4),
        "threshold_averaged_MC_balanced": round(float(((recall + reject) / 2).mean()), 4),
        "R1_crossover_structural_weight": round(float(r1_crossover(recall, reject, tags)), 3),
        "R1_benchmark_actual_weight": round(float(com.truth_struct.mean()), 3),
        "R2_lability": r2,
        "R3_model_free": {
            "definition": f"real destabilizer = |measured| >= {REAL_DESTAB_KCAL} kcal/mol",
            "n_real_destabilizers": n_real, "n_calibration_points": n_cal,
        },
        "gate_fragility": frag,
        "tier_gradient": grad,
        "tradeoff_bootstrap": boot,
    }

    ops.to_csv(out / "COMAVI_audit_operating_points.csv", index=False)
    com2[["variant", "system", "expected_mech_class", "truth_struct", "labile",
          "max_abs_ddg"] + gc].to_csv(out / "COMAVI_audit_lability.csv", index=False)
    inv.to_csv(out / "COMAVI_audit_calibration_inversion.csv", index=False)
    mf.to_csv(out / "COMAVI_audit_model_free_recall.csv", index=False)
    sweep.to_csv(out / "COMAVI_audit_gate_sweep.csv", index=False)
    loso.to_csv(out / "COMAVI_audit_tier_gradient_loso.csv", index=False)
    # ------------------------------------------------------------------
    # A4. Tier ablation: does the tier change any mechanism GRADE?
    #
    # The mechanism classifier reads 'comavi_tier' to decide contact-driven /
    # burial-driven labels. Blanking it and reclassifying isolates the tier's
    # contribution. Reported as a null result in section 3.7.
    # ------------------------------------------------------------------
    abl = com.copy()
    abl_rows = []
    for _, r in abl.iterrows():
        r_no = r.copy()
        r_no["comavi_tier"] = np.nan
        m_with = ac.classify_mechanism_at(r, partners, 2.5)
        m_without = ac.classify_mechanism_at(r_no, partners, 2.5)
        g_with = ac.grade_mechanism_consistency(
            r, m_with, r["expected_mech_class"], ac.classify_axis_status(r))[0]
        g_without = ac.grade_mechanism_consistency(
            r_no, m_without, r["expected_mech_class"], ac.classify_axis_status(r_no))[0]
        abl_rows.append({
            "variant": r["variant"], "system": r["system"],
            "mech_with_tier": m_with, "mech_without_tier": m_without,
            "label_changed": m_with != m_without,
            "grade_with_tier": GMAP.get(g_with), "grade_without_tier": GMAP.get(g_without),
        })
    ablation = pd.DataFrame(abl_rows)
    ablation["grade_changed"] = (
        ablation.grade_with_tier.fillna(-1) != ablation.grade_without_tier.fillna(-1))
    mc_with = ablation.grade_with_tier.mean()
    mc_without = ablation.grade_without_tier.mean()
    summary["tier_ablation"] = {
        "n_graded": int(ablation.grade_with_tier.notna().sum()),
        "mc_with_tier": round(float(mc_with), 4),
        "mc_without_tier": round(float(mc_without), 4),
        "delta": round(float(mc_with - mc_without), 6),
        "n_labels_changed": int(ablation.label_changed.sum()),
        "n_grades_changed": int(ablation.grade_changed.sum()),
        "note": "Tier changes mechanism VOCABULARY, not mechanism ACCURACY.",
    }
    ablation.to_csv(out / "COMAVI_audit_tier_ablation.csv", index=False)

    # ------------------------------------------------------------------
    # A5. Cluster bootstrap on the recall/rejection tradeoff.
    #
    # Resamples whole systems (the paper's primary convention). Confirms the
    # MAGNITUDE of the tradeoff; the SHAPE is monotone by construction and is
    # therefore not evidence (see R2).
    # ------------------------------------------------------------------
    rng = np.random.default_rng(SEED)
    gcols = [f"g_{t}" for t in tags]
    systems = com.system.unique()
    br, bz = [], []
    for _ in range(N_BOOT):
        pick = rng.choice(systems, len(systems), replace=True)
        d = pd.concat([com[com.system == s] for s in pick])
        s_, z_ = d[d.truth_struct], d[~d.truth_struct]
        br.append([s_[c].mean() for c in gcols] if len(s_) else [np.nan] * len(gcols))
        bz.append([z_[c].mean() for c in gcols] if len(z_) else [np.nan] * len(gcols))
    br, bz = np.array(br, float), np.array(bz, float)
    pct = lambda a, i: [round(float(v), 3) for v in np.nanpercentile(a[:, i], [2.5, 97.5])]
    summary["cluster_bootstrap"] = {
        "n_boot": N_BOOT, "seed": SEED, "resample_unit": "system",
        "recall_t10": round(float(recall[0]), 3), "recall_t10_ci": pct(br, 0),
        "recall_tSAP": round(float(recall[-1]), 3), "recall_tSAP_ci": pct(br, -1),
        "rejection_t10": round(float(reject[0]), 3), "rejection_t10_ci": pct(bz, 0),
        "rejection_tSAP": round(float(reject[-1]), 3), "rejection_tSAP_ci": pct(bz, -1),
        "p_tradeoff": round(float(np.nanmean(
            (br[:, 0] > br[:, -1]) & (bz[:, -1] > bz[:, 0]))), 4),
    }

    # ------------------------------------------------------------------
    # A6. Multiplicity correction across the ten gate-sweep tests.
    # ------------------------------------------------------------------
    from statsmodels.stats.multitest import multipletests
    gate_ps = list(sweep.firing_fisher_p) + list(sweep.silent_fisher_p)
    _, padj, _, _ = multipletests(gate_ps, alpha=0.05, method="fdr_bh")
    sweep["firing_fisher_p_bh"] = np.round(padj[:len(sweep)], 4)
    sweep["silent_fisher_p_bh"] = np.round(padj[len(sweep):], 4)
    summary["gate_fragility"]["bh_firing_adj"] = dict(
        zip(sweep.threshold, [round(float(v), 4) for v in padj[:len(sweep)]]))
    summary["gate_fragility"]["bh_silent_adj"] = dict(
        zip(sweep.threshold, [round(float(v), 4) for v in padj[len(sweep):]]))
    summary["gate_fragility"]["n_surviving_bh_05"] = int((padj < 0.05).sum())
    sweep.to_csv(out / "COMAVI_audit_gate_sweep.csv", index=False)

    with open(out / "COMAVI_audit_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)

    print(ops.to_string(index=False))
    print()
    print(mf.to_string(index=False))
    print()
    print(sweep.to_string(index=False))
    print()
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
