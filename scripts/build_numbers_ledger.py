#!/usr/bin/env python3
"""Regenerate the DERIVED numeric keys of reference_outputs/COMAVI_numbers_ledger.json.

That file is read by audit_readme_claims.py, audit_evidence_claims.py and
verify_stress_statistics.py, but nothing in the repository wrote it — it was
hand-maintained, and after the v7.7 evidence-ledger correction it still carried
the pre-correction per-axis tallies (monomer 21/28, fold 20/26, binding 24/32,
tier 34/47; total 99/133; SA 0.7444; MC 0.7193). audit_readme_claims.py detects
that staleness and fails, which in turn prevented verify_audit_binding.py from
instrumenting the README gate at all — so one hand-maintained file disabled two
gates downstream of it.

Only the keys derivable from the canonical are rewritten. Curated keys (labels,
scope_firing_by_class, scope_note, precision_fire_*) are preserved untouched:
they encode decisions, not measurements, and this script has no authority over
them.

Convention: ALL-ROW (every benchmark row), matching the derivation in
audit_readme_claims.py. The per-axis gate is delegated to
apply_concordance_v5.structural_agreement_by_axis rather than re-implemented.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
LEDGER = REPO / "reference_outputs" / "COMAVI_numbers_ledger.json"
STRESS = REPO / "reference_outputs" / "stress_tests" / "comavi_stress_tests.csv"
STATS_RECORD = (REPO / "reference_outputs" / "stress_tests"
                / "COMAVI_stress_verified_statistics.json")
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as concordance  # noqa: E402

AXES = ("monomer", "fold", "binding", "tier")
GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
PARTNER_JUNK = ("_ci95_", "_distinguishable_")
REFERENCE_T = 2.5


def derive(canonical: pd.DataFrame) -> dict:
    partners = [p for p in concordance.discover_partners(canonical)
                if not any(j in p for j in PARTNER_JUNK)]
    totals = {axis: [0, 0] for axis in AXES}
    for _, row in canonical.iterrows():
        per_axis = concordance.structural_agreement_by_axis(
            row, partners, REFERENCE_T, REFERENCE_T, REFERENCE_T)
        for axis, (num, den) in per_axis.items():
            totals[axis][0] += int(num)
            totals[axis][1] += int(den)

    sa_n = sum(v[0] for v in totals.values())
    sa_d = sum(v[1] for v in totals.values())
    mech = canonical["mech_consistency_t25"].map(GRADE_MAP).dropna()
    silent = canonical["expected_mech_class"].eq("structurally_silent")
    mc = round(float(mech.sum()) / len(mech), 4)

    return {
        "SA_by_axis": {axis: [totals[axis][0], totals[axis][1],
                              round(totals[axis][0] / totals[axis][1], 4)]
                       for axis in sorted(AXES)},
        "SA_total": [sa_n, sa_d],
        "SA": round(sa_n / sa_d, 4),
        "MC": mc,
        "MC_n": int(len(mech)),
        "N_variants": int(len(canonical)),
        "silent_n": int(silent.sum()),
        "canonical_sha256": hashlib.sha256(CANON.read_bytes()).hexdigest(),
    }


def derive_stress(current: dict) -> dict:
    """Stress-layer keys, read from verification/stress_tests.py output.

    Kept separate because it needs that run to have happened; when the output
    is absent these keys are left alone rather than guessed.
    """
    if not STRESS.is_file():
        return {}
    s = pd.read_csv(STRESS).set_index("test")
    _stats = (json.loads(STATS_RECORD.read_text())
              if STATS_RECORD.is_file() else None)

    def interval(row):
        return [float(x) for x in str(s.at[row, "null_mean"]).strip("[]").split(",")]

    def rng(row):
        lo, hi = str(s.at[row, "null_mean"]).split("-")
        return [float(lo), float(hi)]

    out = {
        "stress_permutation": {
            "MC_null_mean": round(float(s.at["permutation null (MC)", "null_mean"]), 4),
            "MC_observed": round(float(s.at["permutation null (MC)", "observed"]), 4),
            "p": float(s.at["permutation null (MC)", "p_value"])},
        "stress_LOSO_MC_range": rng("leave-one-system-out (MC range)"),
        "stress_noise_sd": round(float(s.at["replicate noise (MC)", "null_sd"]), 4),
        # 4-decimal form, taken from the verified statistics record so the
        # ledger and verify_stress_statistics.py cannot disagree on rounding
        # (the 3-decimal reportable form in the summary CSV does not match).
        **({"MC_cluster_ci95": _stats["system_cluster_ci95"]["mc"]["rounded_4"],
            "SA_cluster_ci95": _stats["system_cluster_ci95"]["sa"]["rounded_4"]}
           if _stats else {}),
    }
    # The tier-ablated score is only equal to MC when the ablation changes no
    # grade. That is recorded in the ledger, so it is CHECKED, never assumed.
    if current.get("grades_changed_by_ablation") == 0:
        out["MC_tier_ablated"] = current["MC"]
    return {k: v for k, v in out.items() if v is not None}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if the committed ledger is stale")
    args = ap.parse_args()

    canonical = pd.read_csv(CANON, low_memory=False)
    fresh = derive(canonical)
    fresh.update(derive_stress({**json.loads(LEDGER.read_text()), **fresh}))
    current = json.loads(LEDGER.read_text())

    drift = {k: (current.get(k), v) for k, v in fresh.items() if current.get(k) != v}
    if args.check:
        if drift:
            print(f"FAIL: {len(drift)} derived key(s) stale in {LEDGER.name}")
            for k, (was, now) in drift.items():
                print(f"  {k}: committed {was} -> derived {now}")
            return 1
        print(f"OK: all {len(fresh)} derived keys in {LEDGER.name} match the canonical")
        return 0

    if not drift:
        print("nothing to update; derived keys already current")
        return 0

    for k, (was, now) in drift.items():
        print(f"  {k}: {was} -> {now}")
    current.update(fresh)
    LEDGER.write_text(json.dumps(current, indent=1, sort_keys=True) + "\n")
    print(f"wrote {LEDGER.relative_to(REPO)}  ({len(drift)} key(s) updated; "
          f"curated keys preserved)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
