#!/usr/bin/env python3
"""Sensitivity of the mechanism-pattern score to untested-axis half credit.

WHY THIS EXISTS
---------------
Table 3 grades a call that misses every expected effect as inconsistent. One
class of such misses is a call that asserts an energy change only on axes the
literature never tested: it names no tested effect and contradicts nothing that
was measured. Whether that deserves partial credit is a judgement made after
the results were seen, so it is reported as a sensitivity analysis (S2 Text)
and never enters the primary grading.

This generator regrades every call and every constant strategy with
apply_concordance_v5.grade_mechanism_consistency(..., untested_axis_credit=True)
and reports the score at each threshold and, at 2.5 kcal/mol, the margin over
the strongest constant re-selected in every system-resampled draw. The
resampling (system order, seed, draws) and the statistics are those of
scripts/build_evidence_sensitivity.py, imported rather than copied.

Usage:
    python scripts/build_rubric_sensitivity.py           # write the record
    python scripts/build_rubric_sensitivity.py --check   # verify it is current
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import apply_concordance_v5 as ac  # noqa: E402
import build_evidence_sensitivity as es  # noqa: E402

OUT = REPO / "reference_outputs" / "COMAVI_rubric_sensitivity.json"
TAGS = ["t10", "t15", "t20", "t25", "tSAP"]
THRESHOLDS = {"t10": "1.0", "t15": "1.5", "t20": "2.0", "t25": "2.5", "tSAP": "2.9/2.9/3.5"}


def grade(r, label, axes, credit):
    return ac.grade_mechanism_consistency(r, label, r.expected_mech_class, axes,
                                          untested_axis_credit=credit)[0]


def build():
    g, _led = es.load()
    axes = [ac.classify_axis_status(r) for _, r in g.iterrows()]
    labels = [l for l, a in ac.MECH_LABEL_ASSERTIONS.items() if a[4] != "unevaluable"]
    rows = list(g.iterrows())
    m = np.array([[es.GRADE[grade(r, lab, axes[i], True)] for lab in labels]
                  for i, (_, r) in enumerate(rows)])
    m_primary = np.array([[es.GRADE[grade(r, lab, axes[i], False)] for lab in labels]
                          for i, (_, r) in enumerate(rows)])
    silent = g.silent.to_numpy()
    systems = g.system.to_numpy()
    out = {
        "analysis": ("sensitivity of the mechanism-pattern score to half credit for calls "
                     "that assert an energy change only on axes the literature never tested"),
        "status": "sensitivity analysis only; the primary grading is Table 3 with evidence resolution",
        "files": {"canonical_sha256": es.sha256(es.CANON), "ledger_sha256": es.sha256(es.LEDGER)},
        "constant_labels_graded": len(labels), "draws": es.N_DRAWS, "seed": es.SEED,
        "resampling_unit": "system", "by_threshold": {}, "grade_changes": [],
    }
    best_primary = float(m_primary.mean(axis=0).max())
    best = float(m.mean(axis=0).max())
    for t in TAGS:
        stored = g["mech_consistency_%s" % t].map(es.GRADE).to_numpy(float)
        primary = np.array([es.GRADE[grade(r, r["mech_%s" % t], axes[i], False)]
                            for i, (_, r) in enumerate(rows)])
        assert np.array_equal(stored, primary), "stored %s grades differ from the grader" % t
        sens = np.array([es.GRADE[grade(r, r["mech_%s" % t], axes[i], True)]
                         for i, (_, r) in enumerate(rows)])
        assert (sens >= primary).all(), "the credit can only raise a grade"
        for i in np.flatnonzero(sens != primary):
            r = rows[i][1]
            out["grade_changes"].append({
                "threshold": THRESHOLDS[t], "gene": str(r.gene).upper(), "variant": r.variant,
                "expected_class": r.expected_mech_class, "call": r["mech_%s" % t],
                "primary_grade": r["mech_consistency_%s" % t],
                "sensitivity_grade": grade(r, r["mech_%s" % t], axes[i], True)})
        out["by_threshold"][THRESHOLDS[t]] = {
            "primary_score": round(float(primary.mean()), 4),
            "sensitivity_score": round(float(sens.mean()), 4),
            "sensitivity_weighted_total": float(sens.sum()),
            "margin_over_strongest_constant": round(float(sens.mean()) - best, 4)}
        if t == "t25":
            comavi = sens
    idx = np.arange(len(g))
    st = es.stats(idx, silent, labels, m, comavi)
    uniq = pd.unique(systems)
    rng = np.random.default_rng(es.SEED)
    sc, mg = [], []
    for _ in range(es.N_DRAWS):
        pick = rng.choice(uniq, len(uniq), replace=True)
        bi = np.concatenate([idx[systems == s] for s in pick])
        b = es.stats(bi, silent, labels, m, comavi)
        sc.append(b["score"]); mg.append(b["margin"])
    ci = lambda v: [round(float(x), 4) for x in np.nanpercentile(np.asarray(v, float), [2.5, 97.5])]
    out["reference_2_5"] = {
        "score": round(float(st["score"]), 4), "score_ci95": ci(sc),
        "strongest_constant": round(float(st["best"]), 4),
        "strongest_constant_labels": st["best_labels"],
        "margin": round(float(st["margin"]), 4), "margin_ci95": ci(mg),
        "margin_p_at_or_below_zero": round(float(np.mean(np.asarray(mg) <= 0)), 4)}
    changed = [labels[j] for j in range(len(labels))
               if not np.isclose(m[:, j].mean(), m_primary[:, j].mean())]
    out["constants"] = {
        "strongest_primary": round(best_primary, 4), "strongest_sensitivity": round(best, 4),
        "labels_whose_score_changes": changed}
    assert best <= float(comavi.mean()), "a constant strategy reaches COMAVI under the credit"
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    text = json.dumps(build(), indent=2, ensure_ascii=False) + "\n"
    if args.check:
        if not OUT.exists() or OUT.read_text() != text:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: %s is current" % OUT.name)
        return 0
    OUT.write_text(text)
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
