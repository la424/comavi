#!/usr/bin/env python3
"""Regenerate reference_outputs/stress_tests/COMAVI_stress_verified_statistics.json.

verification/verify_stress_statistics.py gates this record against freshly
computed values and file hashes, but nothing wrote it — it was hand-maintained,
so after the v7.7 correction it still carried the pre-correction observed values
(MC 0.7193, SA 0.75), population totals (99/133 all-row, 99/132 primary) and
cluster intervals (0.625-0.821), and the gate failed on the canonical hash.

Every field is derived from the committed stress draw archive written by
verification/stress_tests.py, except the primary-convention structural-agreement
tally, which is computed from the canonical by delegating to
apply_concordance_v5.structural_agreement_by_axis.

Run verification/stress_tests.py --out-dir reference_outputs/stress_tests first;
this script only transcribes and hashes what that run produced.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
ROOT = REPO / "reference_outputs" / "stress_tests"
OUT = ROOT / "COMAVI_stress_verified_statistics.json"
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
STRESS_SCRIPT = REPO / "verification" / "stress_tests.py"
SUMMARY = ROOT / "comavi_stress_tests.csv"
DRAWS = ROOT / "comavi_stress_draws.npz"

ANALYSIS = "COMAVI primary-continuity stress tests"
AXES = ("monomer", "fold", "binding", "tier")
GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
PARTNER_JUNK = ("_ci95_", "_distinguishable_")
REFERENCE_T = 2.5


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval(draws: np.ndarray) -> dict:
    lo, hi = (float(np.percentile(draws, q)) for q in (2.5, 97.5))
    return {"full_precision": [lo, hi],
            "reportable_3": [round(lo, 3), round(hi, 3)],
            "rounded_4": [round(lo, 4), round(hi, 4)]}


def tally(canonical: pd.DataFrame, graded_only: bool) -> list[int]:
    """Four-output structural agreement; graded_only selects the convention.

    primary = the 57 mechanism-gradeable rows (89/124);
    all-row = every benchmark row (89/125).
    """
    sys.path.insert(0, str(REPO / "scripts"))
    import apply_concordance_v5 as concordance

    partners = [p for p in concordance.discover_partners(canonical)
                if not any(j in p for j in PARTNER_JUNK)]
    frame = (canonical[canonical["mech_consistency_t25"].isin(GRADE_MAP)]
             if graded_only else canonical)
    num = den = 0
    for _, row in frame.iterrows():
        for _axis, (a, b) in concordance.structural_agreement_by_axis(
                row, partners, REFERENCE_T, REFERENCE_T, REFERENCE_T).items():
            num += int(a)
            den += int(b)
    return [num, den]


def build() -> dict:
    for required in (CANON, STRESS_SCRIPT, SUMMARY, DRAWS):
        if not required.is_file():
            raise SystemExit(f"required input missing: {required}")

    z = np.load(DRAWS, allow_pickle=True)
    canonical = pd.read_csv(CANON, low_memory=False)
    scalar = lambda key: z[key].reshape(-1)[0]  # noqa: E731

    return {
        "analysis": ANALYSIS,
        "draws": {"cluster_seed": int(scalar("seed_cluster")),
                  "system_cluster_bootstrap": int(z["cluster_mc"].size),
                  "variant_bootstrap": int(z["boot_mc"].size)},
        "files": {"canonical_sha256": sha256(CANON),
                  "stress_draws_sha256": sha256(DRAWS),
                  "stress_script_sha256": sha256(STRESS_SCRIPT),
                  "stress_summary_sha256": sha256(SUMMARY)},
        "observed": {"mechanism_consistency": float(scalar("mc_obs")),
                     "structural_agreement": float(scalar("sa_obs"))},
        "population": {
            "all_row_structural_agreement": tally(canonical, graded_only=False),
            "canonical_rows": int(len(canonical)),
            "mechanism_graded_variants": int(scalar("mc_n")),
            "primary_structural_agreement": tally(canonical, graded_only=True),
            "system_labels": [str(s) for s in z["cluster_systems"]],
            "systems": int(scalar("cluster_system_count"))},
        "system_cluster_ci95": {"mc": interval(z["cluster_mc"]),
                                "sa": interval(z["cluster_sa"])},
        "variant_bootstrap_ci95": {"mc": interval(z["boot_mc"]),
                                   "sa": interval(z["boot_sa"])},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if the committed record is stale")
    args = ap.parse_args()

    fresh = build()
    if args.check:
        current = json.loads(OUT.read_text()) if OUT.is_file() else {}
        drift = [k for k, v in fresh.items() if current.get(k) != v]
        if drift:
            print(f"FAIL: {len(drift)} stale section(s) in {OUT.name}: "
                  f"{', '.join(drift)}")
            return 1
        print(f"OK: {OUT.name} matches the committed stress draws")
        return 0

    OUT.write_text(json.dumps(fresh, indent=1, sort_keys=True) + "\n")
    obs = fresh["observed"]
    print(f"wrote {OUT.relative_to(REPO)}")
    print(f"  observed MC {obs['mechanism_consistency']:.4f} | "
          f"SA {obs['structural_agreement']:.4f}")
    print(f"  all-row {fresh['population']['all_row_structural_agreement']} | "
          f"primary {fresh['population']['primary_structural_agreement']}")
    print(f"  cluster MC CI {fresh['system_cluster_ci95']['mc']['reportable_3']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
