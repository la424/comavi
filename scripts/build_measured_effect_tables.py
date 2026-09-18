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
    """[B] agreement / borderline / substantive disagreement, with median |error|.

    Two corrections, both found by auditing the manuscript table against the
    figure that plots the same quantity (figures/src/figure_delta_calibration.py,
    panel c):

    [1] POPULATION. The barnase Glu73 buried-charge cluster is a prespecified
        force-field failure and is excluded here, as the table's own caption and
        the figure both state. The first version of this generator used all 63
        rows, which contradicted the caption it was supposed to produce and
        reported 43/8/12 where the correct answer is 38/7/11.

    [2] CLASSIFICATION. `straddles` and `both_near` are columns of the
        calibration table. Re-deriving them from a window around the threshold
        looks equivalent and is not: the re-derivation put 8 comparisons in the
        borderline class against the stored columns' 7. Read the columns.
    """
    t = REFERENCE_T
    keep = ~cal.get("g73", pd.Series(False, index=cal.index)).fillna(False).astype(bool)
    cal = cal.loc[keep]
    err = (cal["foldx_ddg"] - cal["measured_kcal"]).abs()
    straddles = cal["straddles"].fillna(False).astype(bool)
    both_near = cal["both_near"].fillna(False).astype(bool)

    agree = ~straddles
    borderline = straddles & both_near
    substantive = straddles & ~both_near

    out = []
    for name, mask, definition in [
            ("Agreement", agree, "Measured and predicted calls agree."),
            ("Borderline disagreement", borderline,
             f"Calls straddle {t} kcal/mol but both values sit near it."),
            ("Substantive disagreement", substantive,
             f"Calls straddle {t} kcal/mol and at least one value is far from it.")]:
        ex = cal.loc[mask, "variant"].dropna().unique().tolist()
        out.append({"call_relationship": name, "n": int(mask.sum()),
                    "definition": definition,
                    "median_abs_error": round(float(err[mask].median()), 2)
                                        if mask.any() else None,
                    "examples": ex[:3]})
    total = sum(r["n"] for r in out)
    assert total == len(cal), f"classes sum to {total}, not {len(cal)}"
    # Guard both corrections: the Glu73 exclusion (n) and reading the stored
    # columns rather than re-deriving them (the class split).
    assert len(cal) == 56, f"expected 56 non-Glu73 comparisons, got {len(cal)}"
    got = tuple(r["n"] for r in out)
    assert got == (38, 7, 11), f"class split {got} != figure panel c (38, 7, 11)"
    return {"n_comparisons": len(cal), "reference_threshold": t,
            "population": "63 comparisons less the 7-point barnase Glu73 cluster",
            "classes": out}


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
