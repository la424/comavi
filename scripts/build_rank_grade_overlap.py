#!/usr/bin/env python3
"""Do the variants COMAVI ranks highly also get their mechanism graded correctly?

A reviewer asked whether the structural mechanisms recovered in the top 20 were
the same variants whose mechanism pattern COMAVI grades correctly. The paper
reported the two results side by side and never crossed them, so the question
was unanswered. It matters because the two could be independent: ranking a
variant highly while mis-calling its mechanism points a user at the right
variant with the wrong hypothesis, which is a different failure from ranking
badly.

Both quantities already exist -- the priority score in the prioritization
output, the mechanism grade in the canonical -- so this is a join, not a new
model.

USAGE
    python scripts/build_rank_grade_overlap.py
    python scripts/build_rank_grade_overlap.py --check
"""

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

REPO = Path(__file__).resolve().parent.parent
PER_VARIANT = REPO / "reference_outputs" / "isds_v1" / "ISDS_v1_per_variant.csv"
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "isds_v1" / "ISDS_v1_rank_grade_overlap.json"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
TOP_K = 20


def build():
    per = pd.read_csv(PER_VARIANT)
    canonical = pd.read_csv(CANON, low_memory=False)

    frame = per.merge(canonical[["system", "variant", "mech_consistency_t25"]],
                      on=["system", "variant"], how="left")
    assert len(frame) == len(per), "join changed row count"

    frame["rank"] = frame["isds_v1"].rank(ascending=False, method="first")
    frame["in_top_k"] = frame["rank"] <= TOP_K
    frame["grade"] = frame["mech_consistency_t25"].map(GRADE_MAP)
    truth = frame["structural_ground_truth"].astype(bool)

    inside = frame[truth & frame["in_top_k"]]
    outside = frame[truth & ~frame["in_top_k"]]

    def named(sub):
        return [{"system": r.system, "gene": r.gene, "variant": r.variant,
                 "rank": int(r["rank"]), "isds_v1": round(float(r.isds_v1), 4),
                 "expected_mech_class": r.expected_mech_class,
                 "grade_label": r.mech_consistency_t25}
                for _, r in sub.sort_values("rank").iterrows()]

    # Association between priority score and mechanism grade. n is small, so
    # the descriptive contrast carries the claim and the tests are reported for
    # completeness rather than as the finding.
    both = frame[truth & frame["grade"].notna()]
    rho, rho_p = spearmanr(both["isds_v1"], both["grade"])
    if len(inside) and len(outside):
        u, u_p = mannwhitneyu(inside["grade"].dropna(), outside["grade"].dropna(),
                              alternative="greater")
    else:
        u, u_p = float("nan"), float("nan")

    result = {
        "question": ("Are the structural mechanisms recovered in the top %d the same "
                     "variants whose mechanism pattern COMAVI grades correctly?" % TOP_K),
        "population": {
            "n_prioritized": int(len(frame)),
            "n_structural_truth": int(truth.sum()),
            "top_k": TOP_K,
        },
        "recovery": {
            "structural_in_top_k": int(len(inside)),
            "structural_outside_top_k": int(len(outside)),
        },
        "mean_grade": {
            "in_top_k": round(float(inside["grade"].mean()), 4),
            "outside_top_k": round(float(outside["grade"].mean()), 4),
        },
        "grade_breakdown_in_top_k": {
            label: int((inside["mech_consistency_t25"] == label).sum())
            for label in ("consistent", "partial", "inconsistent")
        },
        "grade_breakdown_outside_top_k": {
            label: int((outside["mech_consistency_t25"] == label).sum())
            for label in ("consistent", "partial", "inconsistent")
        },
        "consistent_outside_top_k": named(
            outside[outside["mech_consistency_t25"] == "consistent"]),
        "inconsistent_inside_top_k": named(
            inside[inside["mech_consistency_t25"] == "inconsistent"]),
        "structural_outside_top_k_detail": named(outside),
        "association": {
            "spearman_rho_score_vs_grade": round(float(rho), 4),
            "spearman_p": round(float(rho_p), 4),
            "mannwhitney_u_grade_top_vs_rest": None if pd.isna(u) else float(u),
            "mannwhitney_p_one_sided": None if pd.isna(u_p) else round(float(u_p), 4),
            "note": ("n = %d structural-truth variants with %d outside the top %d; "
                     "the descriptive contrast in mean grade carries the claim."
                     % (int(truth.sum()), int(len(outside)), TOP_K)),
        },
    }
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    result = build()

    if args.check:
        if not OUT.is_file():
            print("FAIL: %s missing" % OUT.name)
            return 1
        stored = json.loads(OUT.read_text())
        if stored != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: rank/grade overlap matches a fresh build")
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1) + "\n")

    rec = result["recovery"]
    mg = result["mean_grade"]
    print("structural mechanisms in top %d: %d of %d"
          % (TOP_K, rec["structural_in_top_k"],
             result["population"]["n_structural_truth"]))
    print("  mean mechanism grade  in top %d: %.3f" % (TOP_K, mg["in_top_k"]))
    print("  mean mechanism grade outside   : %.3f" % mg["outside_top_k"])
    print("  breakdown in top %d: %s" % (TOP_K, result["grade_breakdown_in_top_k"]))
    print("  graded consistent yet outside the top %d: %d"
          % (TOP_K, len(result["consistent_outside_top_k"])))
    print("  graded inconsistent yet inside the top %d: %s"
          % (TOP_K, ", ".join("%s %s" % (v["gene"].upper(), v["variant"])
                              for v in result["inconsistent_inside_top_k"]) or "none"))
    print("  Spearman score vs grade: rho=%.3f p=%.4f"
          % (result["association"]["spearman_rho_score_vs_grade"],
             result["association"]["spearman_p"]))
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
