#!/usr/bin/env python3
"""Does the mechanism-pattern result depend on the weakly evidenced variants?

WHY THIS EXISTS
---------------
Evidence quality in the benchmark differs by direction: a disruptive effect is
usually followed up with a dedicated measurement, a normal one rarely is. The
population (E4) and inferred (E5) classes are therefore almost all neutral, and
most of them belong to non-disruptive variants. The paper says so, which raises
the obvious question of whether the headline depends on those variants.

This answers it on evidence-restricted views of the 56-variant benchmark. Each
view is scored against the strongest structure-blind strategy ON THAT VIEW,
because removing non-disruptive variants changes which constant label is
strongest: comparing a restricted score with the full-benchmark floor would
use the wrong baseline. A constant strategy gives every variant the same
mechanism label; every gradeable label the classifier can emit is tried, and
graded with the shipped grader (apply_concordance_v5).

The mechanism-pattern score, like accuracy, rewards the majority call, and the
views differ widely in their share of non-disruptive variants. Each view
therefore also reports the balanced score, the mean of the mechanism arm and
the non-disruptive arm, on which no constant strategy can exceed 0.5.

Intervals resample whole systems (2,000 draws, seed 8, as in the stress tests),
and the strongest constant is re-selected inside every draw, so the margin
interval covers the selection of the baseline too.

Usage:
    python scripts/build_evidence_sensitivity.py           # write the record
    python scripts/build_evidence_sensitivity.py --check   # verify it is current
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402
import benchmark_population as bp  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
LEDGER = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"
OUT = REPO / "reference_outputs" / "COMAVI_evidence_sensitivity.json"
GRADE = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
SILENT = "structurally_silent"
N_DRAWS = 2000
SEED = 8
STRONG = {"E1", "E2", "E3"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load():
    canon = pd.read_csv(CANON, low_memory=False)
    g = bp.benchmark(canon).reset_index(drop=True)
    assert len(g) == bp.N_BENCHMARK, len(g)
    led = pd.read_csv(LEDGER)
    keys = set(zip(g.system, g.variant))
    led = led[[(s, v) in keys for s, v in zip(led.system, led.variant)]].copy()
    led["cls"] = led.evidence_type.str[:2]
    classes = led.groupby(["system", "variant"]).cls.apply(set).to_dict()
    g["classes"] = [classes.get((s, v), set()) for s, v in zip(g.system, g.variant)]
    g["silent"] = g.expected_mech_class == SILENT
    return g, led


def grade_matrix(g):
    """Grades for COMAVI and for every gradeable constant label, per variant."""
    labels = [l for l, a in ac.MECH_LABEL_ASSERTIONS.items() if a[4] != "unevaluable"]
    axes = [ac.classify_axis_status(r) for _, r in g.iterrows()]
    m = np.zeros((len(g), len(labels)))
    for j, lab in enumerate(labels):
        for i, (_, r) in enumerate(g.iterrows()):
            m[i, j] = GRADE[ac.grade_mechanism_consistency(r, lab, r.expected_mech_class, axes[i])[0]]
    comavi = g.mech_consistency_t25.map(GRADE).to_numpy(float)
    assert not np.isnan(comavi).any()
    return labels, m, comavi


def views(g):
    weak_negative = g.silent & g.classes.map(lambda c: not (c & STRONG))
    e4_only = g.silent & g.classes.map(lambda c: c == {"E4"})
    whole_variant = g.silent & g.classes.map(lambda c: len(c) == 0)
    return {
        "full_benchmark": pd.Series(True, index=g.index),
        "without_population_only_negatives": ~e4_only,
        "without_whole_variant_negatives": ~whole_variant,
        "without_weak_evidence_negatives": ~weak_negative,
        "every_expectation_E1_E3": g.classes.map(lambda c: bool(c) and c <= STRONG),
        "every_expectation_E1_E2": g.classes.map(lambda c: bool(c) and c <= {"E1", "E2"}),
    }


def stats(idx, silent, labels, m, comavi):
    s = silent[idx]
    score = comavi[idx].mean()
    const = m[idx].mean(axis=0)
    best = const.max()
    arm_m = comavi[idx][~s].mean() if (~s).any() else np.nan
    arm_s = comavi[idx][s].mean() if s.any() else np.nan
    c_arm_m = m[idx][~s].mean(axis=0) if (~s).any() else np.full(len(labels), np.nan)
    c_arm_s = m[idx][s].mean(axis=0) if s.any() else np.full(len(labels), np.nan)
    return dict(score=score, best=best, best_labels=[labels[j] for j in np.flatnonzero(const == best)],
                margin=score - best, balanced=np.nanmean([arm_m, arm_s]),
                best_balanced=np.nanmax(np.nanmean(np.vstack([c_arm_m, c_arm_s]), axis=0)),
                arm_mechanism=arm_m, arm_non_disruptive=arm_s)


def build(digits=4):
    """The evidence-sensitivity record. Values are rounded to `digits` decimals
    for the committed record; digits=None returns them at full precision, which
    is what any consumer that re-rounds them (the manuscript, to three
    decimals) must use. Re-rounding a 4-decimal value that ends in 5 can move
    the third decimal (v11: the leave-one-system-out minimum 34/53 = 0.64151
    was stored as 0.6415 and printed as 0.641)."""
    R = (lambda x: float(x)) if digits is None else (lambda x: round(float(x), digits))
    g, led = load()
    labels, m, comavi = grade_matrix(g)
    silent = g.silent.to_numpy()
    systems = g.system.to_numpy()
    out = {"analysis": "evidence-restricted views of the mechanism-pattern score at 2.5 kcal/mol",
           "files": {"canonical_sha256": sha256(CANON), "ledger_sha256": sha256(LEDGER)},
           "constant_labels_graded": len(labels), "draws": N_DRAWS, "seed": SEED,
           "resampling_unit": "system", "views": {}}
    for name, mask in views(g).items():
        idx = np.flatnonzero(mask.to_numpy())
        st = stats(idx, silent, labels, m, comavi)
        # Systems in order of first appearance, as verification/stress_tests.py
        # draws them: with the same seed and order the full-benchmark draws are
        # the stress tests' own, and the score interval reproduces theirs.
        uniq = pd.unique(systems[idx])
        rng = np.random.default_rng(SEED)
        sc, mg, bal = [], [], []
        for _ in range(N_DRAWS):
            pick = rng.choice(uniq, len(uniq), replace=True)
            bi = np.concatenate([idx[systems[idx] == s] for s in pick])
            b = stats(bi, silent, labels, m, comavi)
            sc.append(b["score"]); mg.append(b["margin"]); bal.append(b["balanced"])
        ci = lambda v: [R(x) for x in np.nanpercentile(np.asarray(v, float), [2.5, 97.5])]
        removed = g.loc[~mask, ["system", "gene", "variant"]]
        out["views"][name] = {
            "n": int(len(idx)), "non_disruptive": int(silent[idx].sum()),
            "mechanisms": int((~silent[idx]).sum()), "systems": int(len(uniq)),
            "comavi_score": R(st["score"]),
            "comavi_score_ci95": ci(sc),
            "strongest_constant": R(st["best"]),
            "strongest_constant_labels": st["best_labels"],
            "margin": R(st["margin"]),
            "margin_ci95": ci(mg),
            "margin_p_at_or_below_zero": R(np.mean(np.asarray(mg) <= 0)),
            "comavi_balanced": R(st["balanced"]),
            "comavi_balanced_ci95": ci(bal),
            "best_constant_balanced": R(st["best_balanced"]),
            "arm_mechanism": R(st["arm_mechanism"]),
            "arm_non_disruptive": R(st["arm_non_disruptive"]),
            "removed": [f"{str(r.gene).upper()} {r.variant}" for r in removed.itertuples()],
        }
    # Where the full-benchmark margin comes from. Against the strongest constant
    # on the full benchmark (calling every variant non-disruptive), each system
    # contributes the sum over its variants of COMAVI's grade minus that
    # strategy's grade; the contributions add up to the margin times 56. The
    # leave-one-system-out margin re-selects the strongest constant without the
    # system, as the views above do.
    j0 = labels.index("No structural effect detected")
    allidx = np.arange(len(g))
    by_system = {}
    for s in pd.unique(systems):
        sel = systems == s
        keep = allidx[~sel]
        lo = stats(keep, silent, labels, m, comavi)
        by_system[s] = {
            "variants": int(sel.sum()), "non_disruptive": int(silent[sel].sum()),
            "contribution_vs_calling_nothing": R((comavi[sel] - m[sel, j0]).sum()),
            "margin_without_system": R(lo["margin"]),
            "strongest_constant_without_system": lo["best_labels"][0],
        }
    out["by_system"] = dict(sorted(by_system.items(),
                                   key=lambda kv: kv[1]["contribution_vs_calling_nothing"]))
    total = sum(v["contribution_vs_calling_nothing"] for v in by_system.values())
    full = out["views"]["full_benchmark"]
    assert abs(total / len(g) - full["margin"]) < 1e-3, (total, full["margin"])
    assert full["n"] == 56 and full["non_disruptive"] == 27, full
    stress = json.loads((REPO / "reference_outputs" / "stress_tests"
                         / "COMAVI_stress_verified_statistics.json").read_text())
    ref = [round(float(x), 4) for x in stress["system_cluster_ci95"]["mc"]["full_precision"]]
    out["cross_check"] = {"stress_tests_system_cluster_ci95": ref,
                          "full_benchmark_score_ci95": full["comavi_score_ci95"]}
    assert max(abs(a - b) for a, b in zip(ref, full["comavi_score_ci95"])) < 5e-4, out["cross_check"]
    assert all(v["best_constant_balanced"] <= 0.5 + 1e-9 for v in out["views"].values())
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    fresh = build()
    text = json.dumps(fresh, indent=2, ensure_ascii=False) + "\n"
    if args.check:
        if not OUT.exists() or OUT.read_text() != text:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: %s is current" % OUT.name)
        return 0
    OUT.write_text(text)
    print("%-36s %3s %3s %7s %7s %7s  %-17s %7s" % ("view", "n", "ND", "COMAVI", "const", "margin", "margin CI", "balanced"))
    for k, v in fresh["views"].items():
        print("%-36s %3d %3d %7.3f %7.3f %+7.3f  [%+.3f, %+.3f] %7.3f"
              % (k, v["n"], v["non_disruptive"], v["comavi_score"], v["strongest_constant"],
                 v["margin"], v["margin_ci95"][0], v["margin_ci95"][1], v["comavi_balanced"]))
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
