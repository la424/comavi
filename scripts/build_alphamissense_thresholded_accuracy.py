#!/usr/bin/env python3
"""AlphaMissense thresholded accuracy against the benchmark's clinical labels.

WHY THIS EXISTS
---------------
The manuscript reported "correct on 38/43 confident calls and 40/49 when
ambiguous scores were thresholded at 0.5 and two no-calls were counted as
errors". Those counts were hand-typed: no file in the repository contained them,
and their denominators (43, 49) match none of the paper's other populations,
which made them the last unverified numbers in the draft. They reproduce
exactly, so the claim was correct and only its provenance was missing.

This is a comparator accuracy, not a COMAVI result. It is reported because a
reader comparing the two methods will ask how the pathogenicity predictor does
under its own published decision rule rather than only as a ranking (the ROC AUC
comparison on the 47-variant complete-case population is the primary
cross-method summary).

POPULATION AND DECISION RULE
----------------------------
Population: the 49 benchmark variants carrying an AlphaMissense class. Every one
of them also carries a curated clinical annotation (`phenotype`), which is the
ground truth here -- clinical, NOT the structural ground truth used everywhere
else in the paper. That is why this quantity was unaffected by the v7.6 ledger
correction, which changed structural labels only.

  confident calls   AM class is likely_pathogenic or likely_benign (n = 43).
                    Correct when the call matches the clinical annotation,
                    with gain-of-function counted as pathogenic.
  thresholded       the 4 ambiguous variants are called pathogenic at
                    AM pathogenicity >= 0.5, and the 2 variants with no
                    AlphaMissense score are counted as errors, giving n = 49.

Usage
-----
    PYTHONPATH=. python scripts/build_alphamissense_thresholded_accuracy.py
    PYTHONPATH=. python scripts/build_alphamissense_thresholded_accuracy.py --check
"""
import argparse
import json
import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_alphamissense_thresholded_accuracy.json"

CONFIDENT = ("likely_pathogenic", "likely_benign")
AMBIGUOUS = "ambiguous"
NO_CALL = "unavailable"
AM_THRESHOLD = 0.5


def build(canon):
    a = canon[canon["AM class"].notna()].copy()
    assert a.phenotype.notna().all(), "a clinically labelled variant lacks a phenotype"
    a["clinically_pathogenic"] = a.phenotype.astype(str).str.startswith("pathogenic")

    conf = a[a["AM class"].isin(CONFIDENT)].copy()
    conf["am_calls_pathogenic"] = conf["AM class"] == "likely_pathogenic"
    conf["correct"] = conf.am_calls_pathogenic == conf.clinically_pathogenic

    amb = a[a["AM class"] == AMBIGUOUS].copy()
    amb["am_calls_pathogenic"] = amb["AM pathogenicity"] >= AM_THRESHOLD
    amb["correct"] = amb.am_calls_pathogenic == amb.clinically_pathogenic

    no_call = a[a["AM class"] == NO_CALL]

    def detail(frame):
        return [{"gene": str(r.gene).upper(), "variant": r.variant,
                 "am_class": r["AM class"],
                 "am_pathogenicity": (None if pd.isna(r["AM pathogenicity"])
                                      else round(float(r["AM pathogenicity"]), 3)),
                 "clinical_annotation": r.phenotype}
                for _, r in frame.iterrows()]

    thresholded_correct = int(conf.correct.sum() + amb.correct.sum())
    return {
        "ground_truth": "curated clinical annotation (phenotype); gain-of-function counted as pathogenic",
        "note": ("Clinical ground truth, not the structural ground truth used for "
                 "COMAVI's mechanism metrics. Unaffected by the v7.6 ledger "
                 "correction, which changed structural labels only."),
        "population": {
            "with_am_class": int(len(a)),
            "clinically_pathogenic": int(a.clinically_pathogenic.sum()),
            "clinically_benign": int((~a.clinically_pathogenic).sum()),
            "am_class_counts": {k: int(v) for k, v in a["AM class"].value_counts().items()},
        },
        "confident_calls": {
            "n": int(len(conf)), "correct": int(conf.correct.sum()),
            "cell": "%d/%d" % (int(conf.correct.sum()), len(conf)),
            "accuracy": round(float(conf.correct.mean()), 3),
            "errors": detail(conf[~conf.correct]),
        },
        "thresholded_all": {
            "n": int(len(a)), "correct": thresholded_correct,
            "cell": "%d/%d" % (thresholded_correct, len(a)),
            "accuracy": round(thresholded_correct / len(a), 3),
            "threshold": AM_THRESHOLD,
            "ambiguous_resolved": detail(amb),
            "ambiguous_correct": int(amb.correct.sum()),
            "no_calls_counted_as_errors": detail(no_call),
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    result = build(canon)

    c, t = result["confident_calls"], result["thresholded_all"]
    assert c["n"] + len(t["ambiguous_resolved"]) + len(t["no_calls_counted_as_errors"]) == t["n"]
    # the values the manuscript states
    assert c["cell"] == "38/43", c["cell"]
    assert t["cell"] == "40/49", t["cell"]

    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: AlphaMissense thresholded accuracy matches a fresh build")
        return 0

    p = result["population"]
    print("population: %d variants with an AlphaMissense class "
          "(%d clinically pathogenic, %d benign)"
          % (p["with_am_class"], p["clinically_pathogenic"], p["clinically_benign"]))
    print("  class counts: %s" % p["am_class_counts"])
    print("\nconfident calls   %s = %.3f" % (c["cell"], c["accuracy"]))
    for e in c["errors"]:
        print("    %-6s %-8s %-18s AM %.3f  clinical %s"
              % (e["gene"], e["variant"], e["am_class"], e["am_pathogenicity"],
                 e["clinical_annotation"]))
    print("\nthresholded at %.1f, no-calls as errors   %s = %.3f"
          % (t["threshold"], t["cell"], t["accuracy"]))
    print("  ambiguous resolved correctly: %d of %d"
          % (t["ambiguous_correct"], len(t["ambiguous_resolved"])))
    for e in t["ambiguous_resolved"]:
        print("    %-6s %-8s AM %.3f  clinical %s"
              % (e["gene"], e["variant"], e["am_pathogenicity"], e["clinical_annotation"]))
    print("  no-calls: %s"
          % ", ".join("%s %s" % (e["gene"], e["variant"])
                      for e in t["no_calls_counted_as_errors"]))

    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
