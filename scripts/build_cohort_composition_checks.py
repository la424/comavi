#!/usr/bin/env python3
"""Do the headline results depend on how the cohort was assembled?

WHY THIS EXISTS
---------------
Two features of the benchmark's construction can shape its results, and neither
is visible in the headline numbers.

1. Structure source. Seven systems use experimental complex structures and seven
   use AlphaFold 3 models. If COMAVI performs differently on the two, the pooled
   mechanism-pattern score mixes two regimes, and a reader planning to run COMAVI
   on predicted complexes needs the predicted-structure number.

2. Clinical labels. Every benign variant in the pathogenicity comparison entered
   the benchmark as a curated no-lesion control, so a structural score separates
   benign variants from pathogenic structural mechanisms partly by construction.
   The informative quantity is how the PATHOGENIC class divides between
   structural mechanisms and variants with no lesion on the modeled axes: that
   division is direct evidence that structural disruption and pathogenicity are
   different properties.

Both are reported here from committed records. Nothing is fitted.

DEFINITIONS
-----------
Structure source is read from the canonical structure_source column. AlphaFold 3
server models are named fold_<job>_model_<n>.pdb; every other file is a Protein
Data Bank entry. The classification refers to the complex structure. In 11
systems isolated-subunit stability uses a separate single-chain AlphaFold 3 model
whatever the complex source (configs/benchmark_systems.yaml), so this split
tests the source of the assembled context, not of every axis.

The difference in mechanism-pattern score (experimental minus AlphaFold 3) gets
a 95% interval from resampling whole systems within each source, 2,000 draws,
the convention of the main uncertainty analysis. With seven systems per source,
model quality cannot be separated from system difficulty, and the record says so.

Usage (repo root)
-----------------
    PYTHONPATH=.:scripts python scripts/build_cohort_composition_checks.py
    PYTHONPATH=.:scripts python scripts/build_cohort_composition_checks.py --check
"""
import argparse
import json
import pathlib
import re

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"
CANON = RO / "scored_61var_canonical.csv"
PRIORITY = RO / "isds_v1" / "ISDS_v1_per_variant.csv"
AM_SET = RO / "isds_v1" / "ISDS_v1_alphamissense_common_set.csv"
OUT = RO / "COMAVI_cohort_composition_checks.json"

GRADE = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
SILENT = "structurally_silent"
N_BOOT, SEED = 2000, 20260923          # as scripts/build_threshold_sweep.py
AF3_FILE = re.compile(r"^fold_.+_model_\d+\.pdb$")
PDB_FILE = re.compile(r"^[0-9][A-Za-z0-9]{3}(?:_[A-Za-z0-9]+)*\.pdb$")


def source_of(fname):
    f = str(fname).strip()
    if AF3_FILE.match(f):
        return "AlphaFold 3"
    if PDB_FILE.match(f):
        return "experimental"
    raise ValueError("cannot classify structure file %r" % f)


def auc(y, s):
    """ROC AUC as the Mann-Whitney probability, ties counted one half."""
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    pos, neg = s[y], s[~y]
    if not len(pos) or not len(neg):
        return None
    wins = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return round(float(wins / (len(pos) * len(neg))), 4)


def build():
    can = pd.read_csv(CANON, low_memory=False)
    can["source"] = can.structure_source.map(source_of)
    per_system = can.groupby("system").source.nunique()
    assert (per_system == 1).all(), per_system[per_system != 1]
    sys_source = can.groupby("system").source.first()

    g = can[can.mech_consistency_t25.isin(GRADE)].copy()
    g["score"] = g.mech_consistency_t25.map(GRADE)
    by_source = {}
    for src, sub in g.groupby("source"):
        by_source[src] = {"systems": int(sub.system.nunique()), "variants": int(len(sub)),
                          "score_total": float(sub.score.sum()),
                          "score": round(float(sub.score.mean()), 4)}

    # System resampling within each source.
    pools = {src: [sub.score.to_numpy() for _, sub in part.groupby("system")]
             for src, part in g.groupby("source")}
    rng = np.random.default_rng(SEED)
    diffs = np.empty(N_BOOT)
    for b in range(N_BOOT):
        means = {}
        for src, systems in pools.items():
            pick = rng.integers(0, len(systems), len(systems))
            means[src] = np.concatenate([systems[i] for i in pick]).mean()
        diffs[b] = means["experimental"] - means["AlphaFold 3"]
    lo, hi = np.percentile(diffs, [2.5, 97.5])

    pv = pd.read_csv(PRIORITY)
    pv["source"] = pv.system.map(sys_source)
    assert pv.source.notna().all()
    prio = {}
    for src, sub in pv.groupby("source"):
        y = sub.structural_ground_truth.astype(bool)
        prio[src] = {"variants": int(len(sub)), "structural": int(y.sum()),
                     "roc_auc": auc(y, sub.isds_v1)}

    am = pd.read_csv(AM_SET)
    keyed = can.set_index(["system", "variant"])
    graded = keyed.mech_consistency_t25.isin(GRADE)

    def curated_class(system, variant):
        if not graded.loc[(system, variant)]:
            return "ungraded"
        return "no_lesion" if keyed.loc[(system, variant), "expected_mech_class"] == SILENT \
            else "structural"

    am["curated"] = [curated_class(s, v) for s, v in zip(am.system, am.variant)]
    am["label"] = am.clinical_y.map({1: "pathogenic", 0: "benign"})
    assert am.label.notna().all()
    table = {lab: {c: int(((am.label == lab) & (am.curated == c)).sum())
                   for c in ("structural", "no_lesion", "ungraded")}
             for lab in ("pathogenic", "benign")}
    benign_all_no_lesion = table["benign"]["structural"] == 0 and table["benign"]["ungraded"] == 0

    return {
        "structure_source": {
            "rule": ("complex structure file: fold_<job>_model_<n>.pdb is an AlphaFold 3 "
                     "server model; any other file is a Protein Data Bank entry"),
            "systems": {k: v for k, v in sys_source.sort_index().items()},
            "mechanism_pattern_score_at_2.5": by_source,
            "difference_experimental_minus_af3": {
                "value": round(by_source["experimental"]["score"]
                               - by_source["AlphaFold 3"]["score"], 4),
                "ci95_system_resampling": [round(float(lo), 4), round(float(hi), 4)],
                "draws": N_BOOT, "seed": SEED},
            "prioritization": prio,
            "note": ("Seven systems per source: model quality cannot be separated from "
                     "system difficulty. Isolated-subunit stability uses single-chain "
                     "AlphaFold 3 models in 11 systems whatever the complex source."),
        },
        "clinical_label_by_curated_class": {
            "population": "AlphaMissense comparison set (%d variants)" % len(am),
            "table": table,
            "benign_all_curated_no_lesion": bool(benign_all_no_lesion),
            "note": ("Benign variants entered the benchmark as curated no-lesion "
                     "controls, so pathogenicity AUCs of structural scores partly "
                     "reflect construction. The division of the pathogenic class "
                     "is the informative quantity."),
        },
    }


def same(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) < 1e-9
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return a == b


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="verify the committed record matches a fresh build; write nothing")
    args = ap.parse_args()
    res = build()
    s = res["structure_source"]
    assert sum(v["systems"] for v in s["mechanism_pattern_score_at_2.5"].values()) == 14
    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing -- run without --check" % OUT.name)
            return 1
        if not same(json.loads(OUT.read_text()), res):
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: %s reproduces from the committed records" % OUT.name)
        return 0
    OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    for src, v in s["mechanism_pattern_score_at_2.5"].items():
        print("  %-12s %d systems, %d variants, score %.1f/%d = %.3f | prioritization %s"
              % (src, v["systems"], v["variants"], v["score_total"], v["variants"], v["score"],
                 s["prioritization"][src]))
    d = s["difference_experimental_minus_af3"]
    print("  difference %.3f (95%% CI %.3f to %.3f)" % (d["value"], *d["ci95_system_resampling"]))
    print("  clinical label by curated class:", res["clinical_label_by_curated_class"]["table"])
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
