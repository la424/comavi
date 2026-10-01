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
    """[A] measured destabilizations recovered at each threshold.

    v7.14 correction. The earlier version took absolute values on BOTH sides:
    the population was |measured| >= 1.0 and a comparison counted as recovered
    when |predicted| >= t. That counted sign-inverted predictions as recoveries
    -- all seven barnase Glu73 substitutions (measured destabilizing, predicted
    stabilizing) at 1.0 kcal/mol, and TEM1-BLIP YB50A (measured stabilizing,
    predicted destabilizing) at every threshold -- and it put two measured
    STABILIZING comparisons into a population the manuscript calls "measured
    destabilizing effects". The manuscript's own definition says a recovery
    must get the direction right. The population is now measured >= +1.0 and a
    recovery requires predicted >= +t, i.e. the same sign as the measurement.

    The primary result uses every comparison. The Glu73 cluster, a single-site
    systematic sign inversion excluded from nothing by design, is reported as
    a sensitivity view rather than removed.
    """
    g73 = cal.get("g73", pd.Series(False, index=cal.index)).fillna(False).astype(bool)
    inb = cal["in_benchmark"].astype(str).str.lower().isin(("true", "yes"))
    destab = cal["measured_kcal"] >= DESTABILIZER_MIN

    def block(mask):
        n = int(mask.sum())
        rows = []
        for label, t in THRESHOLDS:
            k = int((mask & (cal["foldx_ddg"] >= t)).sum())
            rows.append({"threshold": label, "recovered": k, "n": n,
                         "fraction": round(k / n, 3), "cell": f"{k}/{n}"})
        return n, rows

    n, rows = block(destab)
    n_ex, rows_ex = block(destab & ~g73)
    excluded_stabilizing = cal.loc[(cal["measured_kcal"] <= -DESTABILIZER_MIN), "variant"].tolist()
    return {"population": {"n": n, "in_benchmark": int((destab & inb).sum()),
                           "external": int((destab & ~inb).sum()),
                           "definition": (f"measured >= +{DESTABILIZER_MIN} kcal/mol; recovered when "
                                          "predicted >= +threshold (same sign)"),
                           "measured_stabilizing_not_in_population": excluded_stabilizing},
            "by_threshold": rows,
            "sensitivity_excluding_glu73": {"n": n_ex, "by_threshold": rows_ex}}


def _classes(cal: pd.DataFrame) -> list:
    t = REFERENCE_T
    err = (cal["foldx_ddg"] - cal["measured_kcal"]).abs()
    straddles = cal["straddles"].fillna(False).astype(bool)
    both_near = cal["both_near"].fillna(False).astype(bool)
    out = []
    for name, mask, definition in [
            ("Agreement", ~straddles, "Measured and predicted calls agree."),
            ("Borderline disagreement", straddles & both_near,
             f"Calls straddle {t} kcal/mol but both values sit near it."),
            ("Substantive disagreement", straddles & ~both_near,
             f"Calls straddle {t} kcal/mol and at least one value is far from it.")]:
        ex = cal.loc[mask, "variant"].dropna().unique().tolist()
        out.append({"call_relationship": name, "n": int(mask.sum()),
                    "definition": definition,
                    "median_abs_error": round(float(err[mask].median()), 2)
                                        if mask.any() else None,
                    "examples": ex[:3]})
    assert sum(r["n"] for r in out) == len(cal)
    return out


def call_relationships(cal: pd.DataFrame) -> dict:
    """[B] agreement / borderline / substantive disagreement, with median |error|.

    Classification reads the stored `straddles` and `both_near` columns;
    re-deriving them from a window looks equivalent and is not (it moved one
    comparison between classes). v7.14: the primary result uses all 63
    comparisons; the Glu73 cluster is a named sensitivity view, not an
    exclusion. The exclusion was introduced with the analysis, so it is not
    described as prespecified anywhere.
    """
    g73 = cal.get("g73", pd.Series(False, index=cal.index)).fillna(False).astype(bool)
    return {"n_comparisons": int(len(cal)), "reference_threshold": REFERENCE_T,
            "population": "all comparisons in COMAVI_delta_calibration_points.csv",
            "classes": _classes(cal),
            "sensitivity_excluding_glu73": {"n_comparisons": int((~g73).sum()),
                                            "classes": _classes(cal.loc[~g73])}}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=pathlib.Path,
                    default=REPO / "reference_outputs")
    ap.add_argument("--print-only", action="store_true")
    ap.add_argument("--check", action="store_true",
                    help="verify the committed record matches a fresh build")
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

    if args.check:
        committed = args.out_dir / "COMAVI_measured_effect_tables.json"
        if not committed.exists():
            print("FAIL: %s missing" % committed.name)
            return 1
        if json.loads(committed.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % committed.name)
            return 1
        print("PASS: measured-effect tables match a fresh build")
        return 0

    if not args.print_only:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        out = args.out_dir / "COMAVI_measured_effect_tables.json"
        out.write_text(json.dumps(result, indent=1) + "\n")
        print(f"\nwrote {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
