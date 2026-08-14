#!/usr/bin/env python
"""Verify the v30 reconciliation numbers against the canonical scored table.

Every number in COMAVI_v30_branch_reconciliation_{summary.json,memo.md} and the
tier-component / four-state CSVs is regenerated here from
reference_outputs/scored_61var_canonical.csv and asserted. Run after any change
to the canonical table, the tier formula, or the concordance gating.

The population definition is load-bearing and easy to get wrong: the tier
population is rows with a non-null tier AND a committed mechanism class
(structural or structurally_silent). Admitting 'structurally_uncommitted' or
'interface_uncommitted_magnitude' rows gives n=49 and silently shifts five
four-state cells. Mirrors scripts/analyze_tier_energy_gating.build_population.
"""
import importlib.util, json, pathlib, sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO = pathlib.Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
STRUCTURAL = ["mixed_structural", "ppi_destab_mechanism", "fold_mechanism"]
SILENT = "structurally_silent"
FIRING = ["concordant_disruption", "ddg_only"]
BOOT_SEED = 20260813
BOOT_N = 20000

spec = importlib.util.spec_from_file_location("ac", REPO / "scripts" / "apply_concordance_v5.py")
ac = importlib.util.module_from_spec(spec); spec.loader.exec_module(ac)

sys.path.insert(0, str(REPO / "scripts"))
from comavi_v7.mechanism import assign_tier, compute_disruption_score

FAILS = []
def check(label, got, want, tol=0.0):
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    print(f"  [{'OK' if ok else 'FAIL'}] {label}: {got} (expected {want})")
    if not ok:
        FAILS.append(label)

def rebuild_tier(row, partners, use_interface_bonus=True):
    """Tier from the SHIPPED scorer, not a reimplementation.

    The ablation zeroes the `multi_*_is_interface` flags, which is exactly the
    interface-bonus input and nothing else: the contact terms read
    `multi_*_inter_contacts` and the burial term reads `multi_*_burial`, both
    untouched. Calling the real scorer means the ablation cannot drift from the
    tier it ablates.
    """
    if not use_interface_bonus:
        row = row.copy()
        for pl in partners:
            col = f"multi_{pl}_is_interface"
            if col in row.index:
                row[col] = False
    return assign_tier(compute_disruption_score(row, partners)["comavi_score"])


def main():
    df = pd.read_csv(CANON)
    partners = ac.discover_partners(df)

    pop = df[df.comavi_tier.notna()
             & df.expected_mech_class.isin(STRUCTURAL + [SILENT])].copy()
    pop["structural_gt"] = pop.expected_mech_class.isin(STRUCTURAL)
    pop["strong_tier"] = pop.comavi_tier.isin(ac.FOOTPRINT_TIERS)

    print("population")
    check("n", len(pop), 47)
    check("structural", int(pop.structural_gt.sum()), 17)
    check("silent", int((~pop.structural_gt).sum()), 30)

    print("\ntier reconstruction (gate: reimplementation must match shipped tiers)")
    rebuilt = pop.apply(lambda r: rebuild_tier(r, partners, True), axis=1)
    check("shipped tiers reproduced", int((rebuilt == pop.comavi_tier).sum()), 47)
    if FAILS:
        print("\nABORT: tier rebuild does not match the shipped tier; "
              "the ablation below would be uninterpretable.")
        return 1

    pop["strong_no_iface"] = pop.apply(
        lambda r: rebuild_tier(r, partners, False), axis=1).isin(ac.FOOTPRINT_TIERS)
    struct = pop[pop.structural_gt]
    def score_of(row, use_bonus=True):
        if not use_bonus:
            row = row.copy()
            for pl in partners:
                col = f"multi_{pl}_is_interface"
                if col in row.index:
                    row[col] = False
        return compute_disruption_score(row, partners)["comavi_score"]

    n_bonus = sum(score_of(r, True) != score_of(r, False) for _, r in struct.iterrows())
    check("structural variants receiving interface bonus", int(n_bonus), 15)

    def screen(col):
        a = int((pop[col] & pop.structural_gt).sum())
        b = int((~pop[col] & pop.structural_gt).sum())
        c = int((pop[col] & ~pop.structural_gt).sum())
        d = int((~pop[col] & ~pop.structural_gt).sum())
        return a, b, c, d, a / (a + b), d / (c + d), fisher_exact([[a, c], [b, d]])[1]

    # Fisher p compared to the full double, not a rounded literal: the manuscript
    # quotes 6.9e-05 and 0.0065, and rounding here before comparing would make the
    # check pass for any nearby value.
    for lab, col, sens, spec, fp in [("full", "strong_tier", 1.0, 0.5667, 6.851195942185343e-05),
                                     ("ablated", "strong_no_iface", 0.8235, 0.6, 0.006546276363860495)]:
        a, b, c, d, s, sp, p = screen(col)
        print(f"\ntier screen ({lab})")
        check(f"{lab} sensitivity", round(s, 4), sens)
        check(f"{lab} specificity", round(sp, 4), spec)
        check(f"{lab} Fisher p", p, fp, tol=abs(fp) * 1e-9)
        check(f"{lab} Fisher p as quoted", f"{p:.2g}", f"{fp:.2g}")

    # Zeroing the interface flag demotes four variants, not three. Three are
    # structural (the sensitivity cost); the fourth, D96V, is structurally
    # silent and is the ENTIRE source of the ablated screen's specificity gain
    # (17/30 -> 18/30). Reporting only the structural three understates what the
    # ablation does.
    dem = pop[pop.strong_tier & ~pop.strong_no_iface]
    check("structural variants demoted", sorted(dem.loc[dem.structural_gt, "variant"]),
          ["E542K", "E545K", "W37Y"])
    check("silent variants demoted", sorted(dem.loc[~dem.structural_gt, "variant"]),
          ["D96V"])

    print("\nsensitivity-equality identity (weak-tier structural cell must be empty)")
    check("weak-tier structural variants", int((~pop.strong_tier & pop.structural_gt).sum()), 0)

    print("\nfour-state sweep")
    expect = {
        "t10":  (16, 9, 0, 8, 1, 4, 0, 9),
        "t15":  (15, 6, 0, 6, 2, 7, 0, 11),
        "t20":  (15, 4, 0, 5, 2, 9, 0, 12),
        "t25":  (15, 3, 0, 3, 2, 10, 0, 14),
        "tSAP": (11, 3, 0, 2, 6, 10, 0, 15),
    }
    for tag, _ in ac.THRESHOLD_SPECS:
        f = pop[f"p1_ddg_concordance_{tag}"].isin(FIRING)
        g, s = pop.structural_gt, pop.strong_tier
        got = (int((f & s & g).sum()), int((f & s & ~g).sum()),
               int((f & ~s & g).sum()), int((f & ~s & ~g).sum()),
               int((~f & s & g).sum()), int((~f & s & ~g).sum()),
               int((~f & ~s & g).sum()), int((~f & ~s & ~g).sum()))
        check(f"{tag} cells", got, expect[tag])

    print("\nspecificity-gain cluster bootstrap (resample systems)")
    sil = pop[~pop.structural_gt]
    fires = sil[f"p1_ddg_concordance_t25"].isin(FIRING)
    point = float((~(fires & sil.strong_tier)).mean() - (~fires).mean())
    check("point specificity gain at t=2.5", round(point, 4), 0.1)

    groups = [g.index.values for _, g in pop.groupby("system")]
    keep = {i: (bool(fires[i]), bool(sil.strong_tier[i])) for i in sil.index}
    rng = np.random.default_rng(BOOT_SEED)
    gains = []
    for _ in range(BOOT_N):
        n = f_ = s_ = 0
        for j in rng.integers(0, len(groups), len(groups)):
            for i in groups[j]:
                if i in keep:
                    fi, si = keep[i]; n += 1
                    f_ += fi; s_ += fi and si
        if n:
            gains.append((n - s_) / n - (n - f_) / n)
    gains = np.array(gains)
    check("bootstrap median", round(float(np.median(gains)), 5), 0.09091)
    check("bootstrap CI low", round(float(np.percentile(gains, 2.5)), 5), 0.0)
    check("bootstrap CI high", round(float(np.percentile(gains, 97.5)), 5), 0.24242)
    check("fraction zero-or-negative gain", round(float((gains <= 0).mean()), 5), 0.11355)
    # These 5-dp values are an exact-reproduction check on THIS script's pinned
    # seed and draw count -- they are not reportable precision and must never be
    # quoted in prose. verify_tier_construction.py owns the reportable figures
    # (2 dp on the interval, whole percent on the fraction); assert the two
    # independent estimators still agree there, so neither can drift unnoticed.
    _sg = json.loads((REPO / "reference_outputs" / "COMAVI_tier_construction.json")
                     .read_text())["specificity_gain"]
    check("agrees with canonical CI high (2 dp)",
          round(float(np.percentile(gains, 97.5)), 2),
          _sg["reference_t25_cluster_bootstrap_ci"][1])
    check("agrees with canonical no-gain pct",
          round(float((gains <= 0).mean()) * 100),
          _sg["fraction_resamples_no_gain_pct_reportable"])

    print("\nheadline agreement decomposition")
    axis_n, axis_d = {}, {}
    for _, r in df.iterrows():
        for k, v in (ac.structural_agreement_by_axis(r, partners, 2.5, 2.5, 2.5) or {}).items():
            if v is None:
                continue
            axis_n[k] = axis_n.get(k, 0) + v[0]
            axis_d[k] = axis_d.get(k, 0) + v[1]
    for k, want in [("tier", (34, 47)), ("monomer", (21, 27)),
                    ("fold", (20, 25)), ("binding", (24, 32))]:
        check(f"{k} axis", (axis_n[k], axis_d[k]), want)
    check("all-rows four-output", (sum(axis_n.values()), sum(axis_d.values())), (99, 131))

    unobs = set(ac.unobservable_variants())
    graded = df[~df.variant.isin(unobs)]
    gn = gd = 0
    for _, r in graded.iterrows():
        got = ac.compute_structural_agreement(r, partners, 2.5, 2.5, 2.5)
        if got:
            gn += got[0]; gd += got[1]
    check("primary four-output (unobservable excluded)", (gn, gd), (99, 130))

    gmap = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
    mc = graded.mech_consistency_t25.map(lambda v: gmap.get(str(v).lower())).dropna()
    check("mechanism consistency", (round(float(mc.sum()), 1), len(mc)), (41.0, 57))

    print("\n" + "=" * 62)
    if FAILS:
        print(f"[FAIL] {len(FAILS)} check(s) failed: {FAILS}")
        return 1
    print("[PASS] every v30 reconciliation number regenerates from the canonical table.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
