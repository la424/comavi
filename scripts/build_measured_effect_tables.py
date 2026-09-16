#!/usr/bin/env python3
"""Generator for the two measured-energy quantities the manuscript reports.

Both were hand-maintained until v7.7 and had no generator, which is why they
drifted out of step with the data at different times and for different reasons.

  [A] "Measured effects recovered" (Tables 1, 2, S5) — at each FoldX calling
      threshold, how many of the measured destabilizers COMAVI would still call.
      Population: comparisons in COMAVI_delta_calibration_points.csv with
      |measured| >= 1.0 kcal/mol (n = 44; 13 in-benchmark, 31 external SKEMPI).
      NOTE: this is COMPARISON-level and derives from measured-vs-predicted ddG,
      neither of which the v7.7 ledger correction touched. It is therefore
      UNCHANGED by that correction — an earlier audit in this project wrongly
      flagged it as moving, having matched the "44" to the 44-row AlphaFold-
      annotated table (comavi_v7_concordance_annotated.csv) instead of to this one.

  [B] Call-relationship classes (Table S6) — agreement / borderline disagreement
      / substantive disagreement over ALL comparisons, with median absolute
      error per class. The manuscript's 38/7/11 sums to 56 against the present
      63 comparisons; the deficit is 7 and the GdmCl-unfolding arm now carries
      10 rows, so that table predates the monomer-fold cohort expansion rather
      than the v7.7 correction.

Neither quantity depends on expected_mech_class, structural_ground_truth or
axis_signature, so neither is affected by a regrade.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
CALIB = REPO / "reference_outputs" / "COMAVI_delta_calibration_points.csv"

# Thresholds as the manuscript's operating-point tables report them.
THRESHOLDS = [("1.0", 1.0), ("1.5", 1.5), ("2.0", 2.0), ("2.5", 2.5), ("tSAP", 2.9)]
REFERENCE_T = 2.5          # the reference calling threshold
BORDERLINE_WINDOW = 1.0    # "within 1.0 kcal/mol of the reference threshold"
DESTABILIZER_MIN = 1.0     # |measured| floor defining the recovered population


def load() -> pd.DataFrame:
    cal = pd.read_csv(CALIB)
    missing = {"measured_kcal", "foldx_ddg", "in_benchmark"} - set(cal.columns)
    if missing:
        raise SystemExit(f"{CALIB.name} missing columns: {sorted(missing)}")
    return cal


def recovered(cal: pd.DataFrame) -> dict:
    """[A] measured effects recovered at each threshold."""
    d = cal["measured_kcal"].abs() >= DESTABILIZER_MIN
    inb = cal["in_benchmark"].astype(str).str.lower().isin(("true", "yes"))
    n = int(d.sum())
    rows = []
    for label, t in THRESHOLDS:
        k = int((d & (cal["foldx_ddg"].abs() >= t)).sum())
        rows.append({"threshold": label, "recovered": k, "n": n,
                     "fraction": round(k / n, 3), "cell": f"{k}/{n}"})
    return {"population": {"n": n, "in_benchmark": int((d & inb).sum()),
                           "external": int((d & ~inb).sum()),
                           "definition": f"|measured| >= {DESTABILIZER_MIN} kcal/mol"},
            "by_threshold": rows}


def call_relationships(cal: pd.DataFrame) -> dict:
    """[B] agreement / borderline / substantive disagreement, with median |error|."""
    t, w = REFERENCE_T, BORDERLINE_WINDOW
    mcall = cal["measured_kcal"].abs() >= t
    pcall = cal["foldx_ddg"].abs() >= t
    err = (cal["foldx_ddg"] - cal["measured_kcal"]).abs()
    near = (cal["measured_kcal"].abs().sub(t).abs() <= w) & \
           (cal["foldx_ddg"].abs().sub(t).abs() <= w)

    agree = mcall == pcall
    borderline = ~agree & near
    substantive = ~agree & ~near

    out = []
    for name, mask, definition in [
            ("Agreement", agree, "Measured and predicted calls agree."),
            ("Borderline disagreement", borderline,
             f"Both values lie within {w:.1f} kcal/mol of {t}."),
            ("Substantive disagreement", substantive,
             f"At least one value lies more than {w:.1f} kcal/mol from {t}.")]:
        ex = cal.loc[mask, "variant"].dropna().unique().tolist()
        out.append({"call_relationship": name, "n": int(mask.sum()),
                    "definition": definition,
                    "median_abs_error": round(float(err[mask].median()), 2)
                                        if mask.any() else None,
                    "examples": ex[:3]})
    total = sum(r["n"] for r in out)
    assert total == len(cal), f"classes sum to {total}, not {len(cal)}"
    return {"n_comparisons": len(cal), "reference_threshold": t, "classes": out}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=pathlib.Path,
                    default=REPO / "reference_outputs")
    ap.add_argument("--print-only", action="store_true")
    args = ap.parse_args()

    cal = load()
    result = {"source": CALIB.name,
              "measured_effects_recovered": recovered(cal),
              "call_relationships": call_relationships(cal)}

    rec = result["measured_effects_recovered"]
    print(f"[A] measured effects recovered  (n = {rec['population']['n']}: "
          f"{rec['population']['in_benchmark']} in-benchmark, "
          f"{rec['population']['external']} external)")
    for r in rec["by_threshold"]:
        print(f"      t={r['threshold']:5} {r['cell']:>7}  ({r['fraction']:.3f})")
    print(f"\n[B] call relationships  (n = {result['call_relationships']['n_comparisons']} "
          f"comparisons, reference t = {REFERENCE_T})")
    for c in result["call_relationships"]["classes"]:
        print(f"      {c['call_relationship']:26} n={c['n']:3}  "
              f"median |error| {c['median_abs_error']} kcal/mol  "
              f"e.g. {', '.join(c['examples'][:2])}")

    if not args.print_only:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        out = args.out_dir / "COMAVI_measured_effect_tables.json"
        out.write_text(json.dumps(result, indent=1) + "\n")
        print(f"\nwrote {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
