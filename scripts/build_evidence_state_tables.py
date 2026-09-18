#!/usr/bin/env python3
"""Generator for the two threshold-by-state supplementary tables.

WHY THIS EXISTS
---------------
S10 Table (four evidence states) and S11 Table (structural-mechanism detection
metrics) were hand-maintained and both went stale on the v7.6 ground-truth
correction. They kept the pre-correction population -- 17 structural positives
and 30 negatives -- against the corrected 20 and 27. Neither was caught by the
stale-token sweeps run during the v7.7 reconciliation, because their cells are
paired fractions ("16/9") and rate triples ("0.941 / 0.433 / 0.485") that match
no searchable literal. They were found only by summing the table rows and
noticing the totals were the old population.

The reference-threshold row of the four-state table is also emitted by
analyze_isds_v1.py as `four_state_agreement`; this generator reproduces it and
asserts agreement, so the two cannot diverge.

WHAT IT COMPUTES
----------------
Population: the 47-variant prioritization subset, the same one S4 Table and the
priority-score analyses use -- variants with both a structural-context tier and
a committed endpoint.

Energetic firing is the shipped concordance call (`concordant_disruption` or
`ddg_only`), read per threshold from the canonical. Structural context is
Tier 1-2. The four states cross those two:

    Convergent    energy fires AND Tier 1-2
    Energy only   energy fires, weaker tier
    Context only  Tier 1-2, energy does not fire
    Neither       neither

For each threshold the tables report, per state, the count of structural-truth
variants and of curated no-lesion variants; and for two screens -- any energetic
firing, and convergent -- sensitivity, specificity and positive predictive value.

Usage
-----
    PYTHONPATH=. python scripts/build_evidence_state_tables.py
    PYTHONPATH=. python scripts/build_evidence_state_tables.py --check
"""
import argparse
import json
import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
TIER_CALLS = REPO / "data" / "analysis_inputs" / "isds_v1" / "tier_comparator_variant_calls.csv"
SUMMARY = REPO / "reference_outputs" / "isds_v1" / "ISDS_v1_summary.json"
OUT = REPO / "reference_outputs" / "isds_v1" / "ISDS_v1_evidence_state_tables.json"

# Firing definition used by every tier verifier in this repository.
FIRING = ("concordant_disruption", "ddg_only")
THRESHOLDS = [("t10", "1.0 default"), ("t15", "1.5"), ("t20", "2.0"),
              ("t25", "2.5 reference"), ("tSAP", "tSAP")]
STATES = ["Convergent", "Energy only", "Context only", "Neither"]


def load():
    calls = pd.read_csv(TIER_CALLS)
    canon = pd.read_csv(CANON, low_memory=False)
    conc = [c for c in canon.columns if c.startswith("p1_ddg_concordance_")]
    merged = calls.merge(canon[["system", "variant"] + conc],
                         on=["system", "variant"], how="left")
    assert len(merged) == len(calls), "merge changed row count"
    truth = merged["structural_ground_truth"].astype(bool)
    if "tier_number" in merged.columns:
        tier = merged["tier_number"]
    else:
        tier = merged["reported_tier"].str.extract(r"(\d)").astype(int)[0]
    return merged, truth, tier <= 2


def build(merged, truth, tier12):
    four, metrics = [], []
    for tag, label in THRESHOLDS:
        fires = merged["p1_ddg_concordance_" + tag].isin(FIRING)
        masks = [fires & tier12, fires & ~tier12, ~fires & tier12, ~fires & ~tier12]
        row = {"operating_point": label, "threshold": tag}
        for name, mask in zip(STATES, masks):
            row[name] = {"structural": int((mask & truth).sum()),
                         "no_lesion": int((mask & ~truth).sum())}
        assert sum(row[s]["structural"] for s in STATES) == int(truth.sum())
        assert sum(row[s]["no_lesion"] for s in STATES) == int((~truth).sum())
        four.append(row)

        mrow = {"operating_point": label, "threshold": tag}
        for name, pos in [("any_energetic_firing", fires),
                          ("convergent_energy_and_tier", fires & tier12)]:
            tp = int((pos & truth).sum())
            fn = int((~pos & truth).sum())
            fp = int((pos & ~truth).sum())
            tn = int((~pos & ~truth).sum())
            mrow[name] = {"tp": tp, "fn": fn, "fp": fp, "tn": tn,
                          "sensitivity": round(tp / (tp + fn), 3),
                          "specificity": round(tn / (fp + tn), 3),
                          "ppv": round(tp / (tp + fp), 3) if tp + fp else None}
        metrics.append(mrow)

    return {"population": {"n": len(merged),
                           "structural": int(truth.sum()),
                           "no_lesion": int((~truth).sum())},
            "firing_definition": list(FIRING),
            "four_evidence_states": four,
            "detection_metrics": metrics}


def cross_check(result):
    """The reference row must agree with analyze_isds_v1.py's own record."""
    stored = json.loads(SUMMARY.read_text()).get("four_state_agreement")
    if not stored:
        return "four_state_agreement absent from ISDS_v1_summary.json"
    ref = next(r for r in result["four_evidence_states"] if r["threshold"] == "t25")
    mine = ([ref[s]["structural"] for s in STATES], [ref[s]["no_lesion"] for s in STATES])
    theirs = (list(stored["structural"]), list(stored["no_lesion"]))
    if mine != theirs:
        return "reference row %s != analyze_isds_v1 record %s" % (mine, theirs)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="verify the committed record matches a fresh build")
    args = ap.parse_args()

    merged, truth, tier12 = load()
    result = build(merged, truth, tier12)

    bad = cross_check(result)
    if bad:
        print("FAIL: %s" % bad)
        return 1

    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: evidence-state tables match a fresh build "
              "and the reference row agrees with analyze_isds_v1")
        return 0

    p = result["population"]
    print("population: n=%d  structural=%d  no_lesion=%d"
          % (p["n"], p["structural"], p["no_lesion"]))
    print("\nS10 Table -- four evidence states (structural / no-lesion)")
    print("  %-15s %s" % ("operating point", "  ".join("%-12s" % s for s in STATES)))
    for r in result["four_evidence_states"]:
        cells = ["%d/%d" % (r[s]["structural"], r[s]["no_lesion"]) for s in STATES]
        print("  %-15s %s" % (r["operating_point"], "  ".join("%-12s" % c for c in cells)))
    print("\nS11 Table -- detection metrics (sensitivity / specificity / PPV)")
    for r in result["detection_metrics"]:
        a, c = r["any_energetic_firing"], r["convergent_energy_and_tier"]
        print("  %-15s any %.3f / %.3f / %.3f    convergent %.3f / %.3f / %.3f"
              % (r["operating_point"], a["sensitivity"], a["specificity"], a["ppv"],
                 c["sensitivity"], c["specificity"], c["ppv"]))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
