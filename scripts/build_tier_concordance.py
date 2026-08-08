#!/usr/bin/env python3
"""Tier x ddG concordance: does the structural-evidence tier agree with the
energy axes, and what does each stratum contain?

Regenerates
  reference_outputs/COMAVI_tier_ddg_concordance.csv
  reference_outputs/COMAVI_tier_ddg_concordance.json

Definitions (all from the shipped pipeline, none re-implemented):
  ddG fires   : any of the three axes destabilizing at t = 2.5 under the shipped
                confidence gates (monomer: ddg_monomer + ddg_monomer_confident;
                per partner: ddg_fold_{p} / ddg_binding_{p} gated on
                ddg_{p}_confident and ddg_binding_{p}_indistinguishable).
  tier strong : comavi_tier in {Tier 1, Tier 2}.
  structural  : the literature-defined expected mechanism class is one of
                ppi_destab_mechanism / mixed_structural / fold_mechanism.
                This is the ground-truth column, not a COMAVI output.

Denominator: the 47 tier-carrying variants with a committed expected mechanism
class. The two 'uncommitted' classes (VHL W117R, troponin R162W) are excluded
because their ground truth does not commit to structural or non-structural, and
the 12 BRCT monomer-fold variants carry no tier (the tier is defined on PPI
systems only).
"""
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact, spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

THRESHOLD = 2.5
STRUCT_CLASSES = {"ppi_destab_mechanism", "mixed_structural", "fold_mechanism"}
STRATA = [
    "both agree (tier 1-2 + ΔΔG fires)",
    "tier only (ΔΔG silent)",
    "ΔΔG only (tier 3-4)",
    "neither",
]


def ddg_axis_fire(row, partners, thr=THRESHOLD):
    """Per-axis destabilizing fire under the shipped confidence gates."""
    mono = fold = bind = False
    mv = row.get("ddg_monomer")
    if pd.notna(mv) and bool(row.get("ddg_monomer_confident", False)):
        mono = float(mv) >= thr
    for p in partners:
        pl = str(p).lower()
        if not bool(row.get(f"ddg_{pl}_confident", False)):
            continue
        fv = row.get(f"ddg_fold_{pl}")
        if pd.notna(fv) and float(fv) >= thr:
            fold = True
        if bool(row.get(f"ddg_binding_{pl}_indistinguishable", False)):
            continue
        bv = row.get(f"ddg_binding_{pl}")
        if pd.notna(bv) and float(bv) >= thr:
            bind = True
    return mono, fold, bind


def stratum_of(tier_strong, fires):
    if tier_strong and fires:
        return STRATA[0]
    if tier_strong and not fires:
        return STRATA[1]
    if fires:
        return STRATA[2]
    return STRATA[3]


def main():
    df = pd.read_csv(REPO / "reference_outputs" / "scored_61var_canonical.csv")
    partners = ac.discover_partners(df)

    rows = []
    for _, r in df.iterrows():
        tier = r.get("comavi_tier")
        cls = r.get("expected_mech_class")
        if not isinstance(tier, str) or not tier.startswith("Tier"):
            continue
        if cls not in STRUCT_CLASSES and cls != "structurally_silent":
            continue
        mono, fold, bind = ddg_axis_fire(r, partners)
        fires = mono or fold or bind
        tier_n = int(tier.split()[-1])
        rows.append(
            dict(
                variant=r["variant"],
                system=r["system"],
                label=r.get("phenotype"),
                tier=tier,
                tier_n=tier_n,
                tier_strong=tier_n <= 2,
                fires=fires,
                structural=cls in STRUCT_CLASSES,
                expected_mech_class=cls,
                stratum=stratum_of(tier_n <= 2, fires),
                grade=r.get(f"mech_consistency_t{int(THRESHOLD * 10)}"),
            )
        )
    per = pd.DataFrame(rows)
    per["pathogenic"] = per.label.astype(str).str.startswith("pathogenic")
    per["gnum"] = per.grade.map({"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0})

    tab = (
        per.groupby("stratum")
        .agg(n=("variant", "size"), structural=("structural", "sum"))
        .reindex(STRATA)
        .reset_index()
    )
    tab["silent"] = tab.n - tab.structural
    tab["precision"] = (tab.structural / tab.n).round(3)
    extra = per.groupby("stratum").agg(pathogenic=("pathogenic", "sum"),
                                       mech_consistency=("gnum", "mean"))
    tab = tab.merge(extra.reset_index(), on="stratum", how="left")
    tab["pathogenic_frac"] = (tab.pathogenic / tab.n).round(3)
    tab["mech_consistency"] = tab.mech_consistency.round(3)

    fired = per[per.fires]
    silent = per[~per.fires]

    def gate(sub):
        a = int(((sub.tier_strong) & (sub.structural)).sum())
        b = int(((sub.tier_strong) & (~sub.structural)).sum())
        c = int(((~sub.tier_strong) & (sub.structural)).sum())
        d = int(((~sub.tier_strong) & (~sub.structural)).sum())
        _, p = fisher_exact([[a, b], [c, d]])
        return dict(strong_structural=a, strong_total=a + b, weak_structural=c,
                    weak_total=c + d, fisher_p=round(float(p), 4))

    per["n_streams"] = per.tier_strong.astype(int) + per.fires.astype(int)
    rho_s, p_s = spearmanr(per.n_streams, per.pathogenic)
    streams = per.groupby("n_streams").agg(n=("variant", "size"),
                                           pathogenic=("pathogenic", "sum"),
                                           structural=("structural", "sum"))

    stats = dict(
        n_variants=int(len(per)),
        threshold=THRESHOLD,
        strata={r.stratum: dict(n=int(r.n), structural=int(r.structural),
                                silent=int(r.silent), precision=float(r.precision),
                                pathogenic=int(r.pathogenic),
                                pathogenic_frac=float(r.pathogenic_frac),
                                mech_consistency=float(r.mech_consistency))
                for r in tab.itertuples()},
        gate_among_ddg_firing=gate(fired),
        gate_among_ddg_silent=gate(silent),
        streams={int(i): dict(n=int(r.n), pathogenic=int(r.pathogenic),
                              structural=int(r.structural))
                 for i, r in streams.iterrows()},
        streams_vs_pathogenic=dict(spearman_rho=round(float(rho_s), 4),
                                   p=round(float(p_s), 5)),
    )

    out = REPO / "reference_outputs"
    tab.to_csv(out / "COMAVI_tier_ddg_concordance.csv", index=False)
    per.to_csv(out / "COMAVI_tier_ddg_concordance_per_variant.csv", index=False)
    with open(out / "COMAVI_tier_ddg_concordance.json", "w") as fh:
        json.dump(stats, fh, indent=1)
    print(tab.to_string(index=False))
    print(json.dumps(stats["gate_among_ddg_firing"], indent=1))
    print(json.dumps(stats["gate_among_ddg_silent"], indent=1))


if __name__ == "__main__":
    main()
