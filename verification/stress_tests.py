#!/usr/bin/env python3
"""Stress tests for the COMAVI canonical benchmark.

Four tests a methods reviewer is entitled to ask for. All four re-derive the
expected mechanism class and per-axis status from the ground-truth columns of
the released canonical table using the pipeline's OWN functions
(`apply_concordance_v5`), then grade them against the STORED predictions. The
predictions are never recomputed; FoldX is not re-run.

  A. Permutation null      — is the headline above chance?
  B. Leave-one-system-out  — does any single gene system carry the result?
  C. Replicate-noise       — do calls survive FoldX's own stochasticity?
  D. Variant bootstrap     — confidence intervals on the headlines.

TWO CONVENTIONS THAT MUST BE HONOURED (both cost real debugging time):

  1. `discover_partners()` run over the wide released table also returns the
     per-partner `*_ci95_*` and `*_distinguishable_*` columns as if they were
     partner chains. They must be filtered out or the structural-agreement
     denominator inflates.

  2. The full 61-row canonical resource is loaded. The score() function applies
     one primary continuity population to both metrics by excluding variants
     returned by unobservable_variants(). This yields 57 whole-variant grades
     and 132 applicable structural outputs. The separate all-row denominator
     audit contains 133 outputs.


The permutation null permutes the ground-truth block AS A UNIT. Shuffling each
ground-truth column independently would fabricate mechanism profiles that never
occur in nature and would give a falsely low null.

Usage:
    python3 verification/stress_tests.py \
        --canonical reference_outputs/scored_61var_canonical.csv \
        --scripts-dir scripts \
        --out-dir verification_output
"""

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Canonical calling thresholds (monomer, fold-in-complex, binding) and the
# mechanism-consistency grade -> numeric map.
T25 = (2.5, 2.5, 2.5)
MCMAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
PARTNER_JUNK = ("_ci95_", "_distinguishable_")
NREP = 5          # FoldX BuildModel replicates behind each mean ddG
SEED_PERM = 20260729
SEED_NOISE = 11
SEED_FRAGILITY = 12
SEED_BOOT = 7
SEED_CLUSTER = 8


def clean_partners(df, ac):
    """Partner chains, with the CI/distinguishability columns filtered out."""
    return [p for p in ac.discover_partners(df)
            if not any(s in p for s in PARTNER_JUNK)]


def interaction_rows(df):
    """Return the complete canonical resource used by the stress tests.

    The legacy function name is retained for compatibility. Population
    selection occurs in score(): both metrics exclude variants returned by
    unobservable_variants(), yielding 57 mechanism grades and 132 applicable
    structural outputs. The all-row 99/133 aggregate is a separate denominator
    audit and is not the stress-test headline.
    """
    return df.copy()


def score(df, ac, partners, tag="t25", thr=T25, recall_mech=False):
    """Re-derive expectations from df's ground truth, grade the predictions.

    By default grades the STORED `mech_<tag>` label — correct for tests that
    perturb ground truth (A, B, D). Test C perturbs the ddG columns instead, so
    it must pass recall_mech=True: otherwise the mechanism label is read from
    the unperturbed stored column and mechanism-consistency cannot move at all
    (the test would report sd = 0 by construction).

    Returns (mc_mean, mc_n, sa_num, sa_den).
    """
    d = df.copy()
    if recall_mech:
        d[f"mech_{tag}"] = d.apply(
            lambda r: ac.classify_mechanism_at(r, partners, thr[0]), axis=1)
    d["expected_mech_class"] = d.apply(ac.derive_expected_mech_class, axis=1)
    axis_status = {i: ac.classify_axis_status(r) for i, r in d.iterrows()}
    # Primary continuity convention: variants whose curated mechanism cannot
    # be represented by the deposited structure are excluded from both metrics.
    # The complete canonical resource contains 99/133 agreeing outputs. R1699Q
    # contributes the sole additional all-row monomer output; excluding it and
    # R1699L yields 99/132 over the same 57 variants used for whole-variant
    # mechanism grading.
    unobs = ac.unobservable_variants()
    grades = []
    for i, r in d.iterrows():
        if r.get("variant") in unobs:
            grades.append(None)
            continue
        g, _, _ = ac.grade_mechanism_consistency(
            r, r.get(f"mech_{tag}"), r["expected_mech_class"], axis_status[i])
        grades.append(g)
    mc = pd.Series(grades).map(MCMAP)
    num = den = 0
    for _, r in d.iterrows():
        if r.get("variant") in unobs:
            continue
        n, dd = ac.compute_structural_agreement(r, partners, *thr)
        num += n
        den += dd
    return mc.mean(), int(mc.notna().sum()), num, den


def sd_pairs(df):
    """(value, sd) column pairs for every perturbable ddG axis."""
    val = [c for c in df.columns
           if (c == "ddg_monomer" or c.startswith("ddg_fold_")
               or c.startswith("ddg_binding_"))
           and not any(s in c for s in ("_ci95_", "_distinguishable_",
                                        "_sd", "_runs", "_confident"))]
    return [(c, c + "_sd") for c in val if c + "_sd" in df.columns]


def permute_block(df, gt_cols, rng):
    """Reassign whole ground-truth profiles across rows."""
    d = df.copy()
    d[gt_cols] = d[gt_cols].values[rng.permutation(len(d))]
    return d


def perturb(df, pairs, rng):
    """Add Gaussian noise at the standard error of the replicate mean."""
    d = df.copy()
    for v, s in pairs:
        sd = pd.to_numeric(d[s], errors="coerce").fillna(0.0).values
        d[v] = pd.to_numeric(d[v], errors="coerce").values + rng.normal(0, sd / np.sqrt(NREP))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--canonical", required=True,
                    help="reference_outputs/scored_61var_canonical.csv")
    ap.add_argument("--scripts-dir", default="scripts")
    ap.add_argument("--out-dir", default="verification_output")
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--n-noise", type=int, default=500)
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()

    sys.path.insert(0, str(Path(args.scripts_dir).resolve()))
    ac = importlib.import_module("apply_concordance_v5")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    canon = pd.read_csv(args.canonical, low_memory=False)
    partners = clean_partners(canon, ac)
    gt_cols = list(ac.GROUND_TRUTH_COLS)
    sub = interaction_rows(canon)
    pairs = sd_pairs(canon)

    print(f"canonical rows {len(canon)} | scored rows {len(sub)} | "
          f"partners {len(partners)} | perturbable axes {len(pairs)}")

    # --- Gate: the scorer must reproduce the released headlines exactly -------
    mc_obs, mc_n, sa_num, sa_den = score(sub, ac, partners)
    sa_obs = sa_num / sa_den
    stored_mc = sub["mech_consistency_t25"].map(MCMAP).mean()
    assert abs(mc_obs - stored_mc) < 1e-9, (
        f"scorer does not reproduce stored grades: {mc_obs} vs {stored_mc}")
    print(f"observed  MC {mc_obs:.4f} (n={mc_n})   SA {sa_num}/{sa_den} = {sa_obs:.4f}   [reproduces stored]")

    # --- A. permutation null -------------------------------------------------
    rng = np.random.default_rng(SEED_PERM)
    mc_null, sa_null = [], []
    for _ in range(args.n_perm):
        m, _, n, d = score(permute_block(sub, gt_cols, rng), ac, partners)
        mc_null.append(m)
        sa_null.append(n / d if d else np.nan)
    mc_null = np.asarray(mc_null)
    sa_null = np.asarray(sa_null, dtype=float)
    p_mc = (np.sum(mc_null >= mc_obs) + 1) / (args.n_perm + 1)
    p_sa = (np.sum(sa_null >= sa_obs) + 1) / (args.n_perm + 1)
    print(f"[A] null MC {mc_null.mean():.4f}+-{mc_null.std():.4f} p={p_mc:.5f} | "
          f"null SA {sa_null.mean():.4f}+-{sa_null.std():.4f} p={p_sa:.5f}")

    # --- B. leave-one-system-out --------------------------------------------
    rows = []
    # v7.3: group on `system`, not `gene`. Two reasons the gene key was wrong:
    #   (1) it is NULL for all 12 BRCA1-BRCT rows, so the fold cohort was
    #       dropped as one unnamed "nan" group rather than as a named system;
    #   (2) it splits two complexes across groups (mlh1_pms2 -> mlh1 + pms2,
    #       pi3k -> pik3ca + pik3r1), so "dropping mlh1" left a pms2 row of the
    #       same complex in place — not a leave-one-SYSTEM-out at all.
    # `system` is populated on every row and partitions the table exactly.
    for gene in sorted(sub["system"].astype(str).unique()):
        keep = sub[sub["system"].astype(str) != gene]
        m, n, N, D = score(keep, ac, partners)
        rows.append(dict(dropped=gene,
                         n_dropped=int((sub["system"].astype(str) == gene).sum()),
                         mc=round(m, 4), mc_n=n, sa=f"{N}/{D}",
                         sa_val=round(N / D, 4)))
    loso = pd.DataFrame(rows).sort_values("mc")
    loso.to_csv(out / "comavi_leave_one_system_out.csv", index=False)
    print(f"[B] MC {loso.mc.min():.4f}-{loso.mc.max():.4f} | "
          f"SA {loso.sa_val.min():.4f}-{loso.sa_val.max():.4f}")

    # --- C. FoldX replicate noise ------------------------------------------
    rng = np.random.default_rng(SEED_NOISE)
    noise = []
    for _ in range(args.n_noise):
        d = perturb(sub, pairs, rng)
        m, _, N, D = score(d, ac, partners, recall_mech=True)
        mt = d.apply(lambda r: ac.classify_mechanism_at(r, partners, 2.5), axis=1)
        noise.append((m, N / D if D else np.nan,
                      int((mt.values != sub["mech_t25"].values).sum())))
    noise = np.asarray(noise, dtype=float)
    print(f"[C] MC {noise[:, 0].mean():.4f}+-{noise[:, 0].std():.4f} | "
          f"SA {noise[:, 1].mean():.4f}+-{noise[:, 1].std():.4f} | "
          f"label flips/draw mean {noise[:, 2].mean():.2f} max {int(noise[:, 2].max())} | "
          f"{100 * np.mean(noise[:, 2] == 0):.0f}% draws label-identical")

    # per-variant fragility
    rng = np.random.default_rng(SEED_FRAGILITY)
    flip = np.zeros(len(sub), dtype=int)
    for _ in range(args.n_noise):
        d = perturb(sub, pairs, rng)
        mt = d.apply(lambda r: ac.classify_mechanism_at(r, partners, 2.5), axis=1)
        flip += (mt.values != sub["mech_t25"].values).astype(int)
    frag = (pd.DataFrame({"system": sub["system"].values,
                          "gene": sub["gene"].values,
                          "variant": sub["variant"].values,
                          "mech_t25": sub["mech_t25"].values,
                          "flip_rate": flip / args.n_noise})
            .sort_values("flip_rate", ascending=False))
    frag.to_csv(out / "comavi_noise_fragility.csv", index=False)
    print(f"    fully stable under noise: {int((frag.flip_rate == 0).sum())}/{len(frag)}")

    # --- D. variant bootstrap ----------------------------------------------
    rng = np.random.default_rng(SEED_BOOT)
    idx = np.arange(len(sub))
    bm, bs = [], []
    for _ in range(args.n_boot):
        b = sub.iloc[rng.choice(idx, len(idx), replace=True)].reset_index(drop=True)
        m, _, N, D = score(b, ac, partners)
        bm.append(m)
        bs.append(N / D if D else np.nan)
    bm = np.asarray(bm)
    bs = np.asarray(bs, dtype=float)
    mc_ci = np.nanpercentile(bm, [2.5, 97.5])
    sa_ci = np.nanpercentile(bs, [2.5, 97.5])
    print(f"[D] MC {mc_obs:.4f} 95% CI [{mc_ci[0]:.3f}, {mc_ci[1]:.3f}] | "
          f"SA {sa_obs:.4f} 95% CI [{sa_ci[0]:.3f}, {sa_ci[1]:.3f}]")

    # --- E. cluster bootstrap (resample whole systems) -----------------------
    # Axes within a system are not independent -- they share a structure, a
    # partner set and a curator. Resampling variants (D) therefore understates
    # the interval. Resampling whole systems is the honest analogue, at the
    # cost of a much coarser resampling unit (14 systems).
    rng = np.random.default_rng(SEED_CLUSTER)
    systems = sub["system"].unique()
    cbm, cbs = [], []
    for _ in range(args.n_boot):
        pick = rng.choice(systems, len(systems), replace=True)
        b = pd.concat([sub[sub.system == s] for s in pick]).reset_index(drop=True)
        m, _, N, D = score(b, ac, partners)
        cbm.append(m)
        cbs.append(N / D if D else np.nan)
    cmc_ci = np.nanpercentile(np.asarray(cbm, dtype=float), [2.5, 97.5])
    csa_ci = np.nanpercentile(np.asarray(cbs, dtype=float), [2.5, 97.5])
    print(f"[E] cluster (n={len(systems)} systems) "
          f"MC 95% CI [{cmc_ci[0]:.3f}, {cmc_ci[1]:.3f}] | "
          f"SA 95% CI [{csa_ci[0]:.3f}, {csa_ci[1]:.3f}]")

    summary = pd.DataFrame([
        dict(test="permutation null (MC)", observed=round(mc_obs, 4),
             null_mean=round(mc_null.mean(), 4), null_sd=round(mc_null.std(), 4),
             p_value=round(p_mc, 5)),
        dict(test="permutation null (SA)", observed=round(sa_obs, 4),
             null_mean=round(sa_null.mean(), 4), null_sd=round(sa_null.std(), 4),
             p_value=round(p_sa, 5)),
        dict(test="leave-one-system-out (MC range)", observed=round(mc_obs, 4),
             null_mean=f"{loso.mc.min():.4f}-{loso.mc.max():.4f}",
             null_sd="", p_value=""),
        dict(test="leave-one-system-out (SA range)", observed=round(sa_obs, 4),
             null_mean=f"{loso.sa_val.min():.4f}-{loso.sa_val.max():.4f}",
             null_sd="", p_value=""),
        dict(test="replicate noise (MC)", observed=round(mc_obs, 4),
             null_mean=round(noise[:, 0].mean(), 4),
             null_sd=round(noise[:, 0].std(), 4), p_value=""),
        dict(test="replicate noise (SA)", observed=round(sa_obs, 4),
             null_mean=round(noise[:, 1].mean(), 4),
             null_sd=round(noise[:, 1].std(), 4), p_value=""),
        dict(test="bootstrap 95% CI (MC)", observed=round(mc_obs, 4),
             null_mean=f"[{mc_ci[0]:.3f}, {mc_ci[1]:.3f}]", null_sd="", p_value=""),
        dict(test="bootstrap 95% CI (SA)", observed=round(sa_obs, 4),
             null_mean=f"[{sa_ci[0]:.3f}, {sa_ci[1]:.3f}]", null_sd="", p_value=""),
        dict(test="cluster bootstrap 95% CI (MC)", observed=round(mc_obs, 4),
             null_mean=f"[{cmc_ci[0]:.3f}, {cmc_ci[1]:.3f}]", null_sd="", p_value=""),
        dict(test="cluster bootstrap 95% CI (SA)", observed=round(sa_obs, 4),
             null_mean=f"[{csa_ci[0]:.3f}, {csa_ci[1]:.3f}]", null_sd="", p_value=""),
    ])
    summary.to_csv(out / "comavi_stress_tests.csv", index=False)

    # Raw draws, so the supplementary figure is reproducible from the repo
    # rather than from an interactive session.
    np.savez(out / "comavi_stress_draws.npz",
             mc_null=mc_null, sa_null=sa_null,
             noise_mc=noise[:, 0], noise_sa=noise[:, 1], noise_flips=noise[:, 2],
             boot_mc=bm, boot_sa=bs,
             cluster_mc=np.asarray(cbm, dtype=float),
             cluster_sa=np.asarray(cbs, dtype=float),
             seed_cluster=np.array([SEED_CLUSTER], dtype=int),
             cluster_system_count=np.array([len(systems)], dtype=int),
             cluster_systems=np.asarray(systems, dtype=str),
             mc_n=np.array([mc_n], dtype=int),
             sa_numerator=np.array([sa_num], dtype=int),
             sa_denominator=np.array([sa_den], dtype=int),
             mc_obs=np.array([mc_obs]), sa_obs=np.array([sa_obs]),
             p_mc=np.array([p_mc]), p_sa=np.array([p_sa]))

    print(f"\nwrote {out}/comavi_stress_tests.csv, "
          f"comavi_leave_one_system_out.csv, comavi_noise_fragility.csv, "
          f"comavi_stress_draws.npz")


if __name__ == "__main__":
    main()
