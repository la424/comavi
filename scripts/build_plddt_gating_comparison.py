#!/usr/bin/env python3
"""What the two pLDDT thresholds actually do, and what "strict" means where.

WHY THIS EXISTS
---------------
The method is described as pLDDT-gated at strict >= 70 and relaxed >= 50, and a
reader naturally expects two sets of headline numbers, one per gate. That is not
how the shipped pipeline works, and three different things in this repository are
called "strict" or "relaxed". This generator states all of them with their
numbers so the distinction can be written down once instead of re-derived.

  [1] THE pLDDT GATE (the one in the method description). Both thresholds are
      active simultaneously in a single pipeline, not as two configurations:
        >= 50  admits a partner interface into the gated set at all
               (comavi_v7/mechanism.py)
        >= 70  upgrades ddG confidence from "medium" to "high"
               (apply_concordance_v5.derive_ddg_confidence)
        < 70   multiplies the tier score by 0.7;  < 50 multiplies it by 0.4
      There is therefore no switch to flip. The equivalent of a strict-gate
      analysis is restricting the cohort to crystal or confident site models,
      which is exactly the set whose site pLDDT clears 70. The manuscript
      already runs that restriction but reports it only qualitatively.

      Note: config.py defines PLDDT_STRICT = 70 and PLDDT_RELAXED = 50 and
      NOTHING READS THEM. The live thresholds are literals at three call sites.
      Treat the constants as dead code, not as the source of truth.

  [2] ddg_<axis>_vote_strict / _vote_relaxed are ENERGY thresholds, not pLDDT:
      strict = |ddG| >= 2.0 kcal/mol, relaxed = |ddG| >= 1.0 kcal/mol
      (apply_concordance_v5.compute_axis_votes).

  [3] concordance_strict / concordance_relaxed are four-way evidence
      strictness across tier, ddG, AlphaMissense and Franklin, where only the
      ddG voter's strictness depends on pLDDT, via the confidence tier in [1]
      (comavi_v7/concordance.py).

Anything that reports "strict" without saying which of the three it means is
ambiguous. This emits all three.

Usage
-----
    PYTHONPATH=. python scripts/build_plddt_gating_comparison.py
    PYTHONPATH=. python scripts/build_plddt_gating_comparison.py --check
"""
import argparse
import json
import pathlib
import re

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_plddt_gating_comparison.json"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
PLDDT_STRICT = 70.0
PLDDT_RELAXED = 50.0
CONFIDENT_SITE = ("crystal", "confident")


def partner_plddt_table(canon):
    """One row per variant x partner interface that has a pLDDT value."""
    cols = [c for c in canon.columns if re.fullmatch(r"multi_\w+_plddt", c)]
    rows = []
    for c in cols:
        partner = c[len("multi_"):-len("_plddt")]
        for i, v in canon[c].items():
            if pd.notna(v):
                rows.append({"system": canon.at[i, "system"],
                             "variant": canon.at[i, "variant"],
                             "partner": partner, "plddt": float(v)})
    return pd.DataFrame(rows)


def build(canon):
    pt = partner_plddt_table(canon)
    strict = pt.plddt >= PLDDT_STRICT
    band = (pt.plddt >= PLDDT_RELAXED) & (pt.plddt < PLDDT_STRICT)
    below = pt.plddt < PLDDT_RELAXED

    graded = canon[canon.mech_consistency_t25.isin(GRADE_MAP)].copy()
    graded["score"] = graded.mech_consistency_t25.map(GRADE_MAP)
    conf = graded.site_plddt_status.isin(CONFIDENT_SITE)

    def summary(frame):
        return {"n": int(len(frame)),
                "total": float(frame.score.sum()),
                "mean": round(float(frame.score.mean()), 4)}

    removed = graded.loc[~conf, ["system", "gene", "variant",
                                 "expected_mech_class", "mech_consistency_t25"]]

    return {
        "plddt_gate": {
            "thresholds": {"strict": PLDDT_STRICT, "relaxed": PLDDT_RELAXED},
            "note": ("Both are active in one pipeline: >=50 admits a partner "
                     "interface, >=70 upgrades ddG confidence to high, and a best "
                     "pLDDT below 70 or 50 multiplies the tier score by 0.7 or 0.4. "
                     "config.PLDDT_STRICT and config.PLDDT_RELAXED are unused."),
            "partner_interfaces": {
                "n": int(len(pt)),
                "at_or_above_strict": int(strict.sum()),
                "in_relaxed_band_50_70": int(band.sum()),
                "below_relaxed": int(below.sum())},
            "relaxed_band_members": pt.loc[band].sort_values("plddt")
                                      .to_dict("records"),
            "site_status_counts": {k: int(v) for k, v in
                                   canon.site_plddt_status.value_counts().items()},
        },
        "cohort_under_each_regime": {
            "relaxed_shipped": summary(graded),
            "strict_confident_site_only": summary(graded[conf]),
            "delta_mean": round(float(graded[conf].score.mean()
                                      - graded.score.mean()), 4),
            "variants_removed_by_strict": removed.to_dict("records"),
        },
        "other_strict_relaxed_senses": {
            "per_axis_energy_votes": {
                "meaning": "|ddG| threshold, NOT pLDDT",
                "strict_kcal_mol": 2.0, "relaxed_kcal_mol": 1.0,
                "counts": {axis: {
                    gate: int(canon["ddg_%s_vote_%s" % (axis, gate)].fillna(False)
                              .astype(bool).sum())
                    for gate in ("strict", "relaxed")}
                    for axis in ("monomer", "fold", "binding")
                    if "ddg_%s_vote_strict" % axis in canon.columns},
                "evaluable": {axis: int(canon["ddg_%s_vote_strict" % axis].notna().sum())
                              for axis in ("monomer", "fold", "binding")
                              if "ddg_%s_vote_strict" % axis in canon.columns},
            },
            "four_way_concordance": {
                "meaning": ("evidence strictness across tier, ddG, AlphaMissense "
                            "and Franklin; only the ddG voter depends on pLDDT"),
                "mean_fraction": {
                    gate: round(float(pd.to_numeric(
                        canon["concordance_%s_full" % gate].astype(str)
                        .str.split("/").str[0], errors="coerce").sum()
                        / pd.to_numeric(
                        canon["concordance_%s_full" % gate].astype(str)
                        .str.split("/").str[1], errors="coerce").sum()), 4)
                    for gate in ("strict", "relaxed")
                    if "concordance_%s_full" % gate in canon.columns},
            },
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    result = build(canon)

    g = result["plddt_gate"]["partner_interfaces"]
    assert g["at_or_above_strict"] + g["in_relaxed_band_50_70"] + g["below_relaxed"] == g["n"]
    c = result["cohort_under_each_regime"]
    assert c["relaxed_shipped"]["n"] == 57, c["relaxed_shipped"]["n"]
    assert abs(c["relaxed_shipped"]["total"] - 39.5) < 1e-9

    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: pLDDT gating comparison matches a fresh build")
        return 0

    print("[1] pLDDT gate -- %d partner interfaces" % g["n"])
    print("      >= 70 (strict)      %2d" % g["at_or_above_strict"])
    print("      50-70 (relaxed only) %2d" % g["in_relaxed_band_50_70"])
    print("      < 50                 %2d" % g["below_relaxed"])
    print("    site status: %s" % result["plddt_gate"]["site_status_counts"])
    print("\n    the %d relaxed-only interfaces:" % g["in_relaxed_band_50_70"])
    for r in result["plddt_gate"]["relaxed_band_members"]:
        print("      %-14s %-8s %-10s %.2f" % (r["system"], r["variant"],
                                               r["partner"], r["plddt"]))
    print("\n    headline under each regime:")
    for k in ("relaxed_shipped", "strict_confident_site_only"):
        s = c[k]
        print("      %-28s %.1f/%d = %.4f" % (k, s["total"], s["n"], s["mean"]))
    print("      delta %+.4f" % c["delta_mean"])
    print("\n[2] per-axis ENERGY votes (|ddG| >= 2.0 strict / >= 1.0 relaxed)")
    ev = result["other_strict_relaxed_senses"]["per_axis_energy_votes"]
    for axis, d in ev["counts"].items():
        print("      %-8s strict %2d   relaxed %2d   of %2d evaluable"
              % (axis, d["strict"], d["relaxed"], ev["evaluable"][axis]))
    print("\n[3] four-way concordance (evidence strictness, pLDDT only via confidence)")
    for gate, v in result["other_strict_relaxed_senses"]["four_way_concordance"]["mean_fraction"].items():
        print("      %-8s pooled fraction %.4f" % (gate, v))

    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
