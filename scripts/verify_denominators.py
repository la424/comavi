"""Regenerate every denominator in the COMAVI paper and check the review's claims.

Written to adjudicate an external review that alleged the BRCT rows were graded
under a different rubric and lack the uncertainty-gating fields. Both claims are
testable against the shipped pipeline; this script tests them.

Usage:
    python3 scripts/verify_denominators.py
Exit 0 if every assertion holds.
"""
import sys
import pathlib
import importlib.util

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
CALIB = REPO / "reference_outputs" / "COMAVI_delta_calibration_points.csv"
GMAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}


def load_ac():
    spec = importlib.util.spec_from_file_location(
        "ac", REPO / "scripts" / "apply_concordance_v5.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ac"] = mod
    spec.loader.exec_module(mod)
    return mod


def main():
    ac = load_ac()
    df = pd.read_csv(CANON)
    partners = [p for p in ac.discover_partners(df) if p]
    brct = df.system == "brca1_brct"
    graded = df["mech_consistency_t25"].isin(GMAP)
    fail = []

    def check(label, got, want):
        ok = got == want
        print(f"  {'ok ' if ok else 'FAIL'}  {label:52} {got}")
        if not ok:
            fail.append(f"{label}: got {got}, expected {want}")

    print("\n[1] Cohort composition")
    check("benchmark rows", len(df), 61)
    check("systems", int(df.system.nunique()), 14)
    check("BRCT rows", int(brct.sum()), 12)
    check("graded variants", int(graded.sum()), 57)
    check("graded interaction", int((graded & ~brct).sum()), 47)
    check("graded BRCT", int((graded & brct).sum()), 10)

    print("\n[2] Mechanism-consistency by population")
    for label, mask, want in [
            ("all gradeable n=57", graded, 0.7193),
            ("interaction n=47", graded & ~brct, 0.7234),
            ("BRCT core/fold n=10", graded & brct, 0.7000)]:
        s = df.loc[mask, "mech_consistency_t25"].map(GMAP)
        check(label, round(float(s.mean()), 4), want)

    print("\n[3] Per-axis decomposition reconciles across cohorts")
    def tally(sub):
        tot = {}
        for _, row in sub.iterrows():
            for k, (n, d) in ac.structural_agreement_by_axis(
                    row, partners, 2.5, 2.5, 2.5).items():
                a, b = tot.get(k, (0, 0))
                tot[k] = (a + n, b + d)
        return tot

    t_all, t_int, t_brct = tally(df), tally(df[~brct]), tally(df[brct])
    for k in t_all:
        check(f"{k}: interaction + BRCT == all",
              (t_int[k][0] + t_brct[k][0], t_int[k][1] + t_brct[k][1]), t_all[k])
    check("aggregate structural agreement",
          (sum(v[0] for v in t_all.values()), sum(v[1] for v in t_all.values())),
          (99, 131))
    check("interaction-only structural agreement",
          (sum(v[0] for v in t_int.values()), sum(v[1] for v in t_int.values())),
          (92, 120))

    print("\n[4] REVIEW CLAIM: BRCT graded under a binary convention")
    # The ungraded-by-construction set is removed upstream of the rubric, by
    # curated role -- not by grade_mechanism_consistency. Compare only the rows
    # the pipeline actually grades; the excluded pair is examined in [8].
    excluded = ac.unobservable_variants()
    mismatch, regrades = [], {}
    for _, row in df[brct].iterrows():
        mech = ac.classify_mechanism_at(row, partners, 2.5)
        grade = ac.grade_mechanism_consistency(
            row, mech[0] if isinstance(mech, tuple) else mech,
            row.get("expected_mech_class"), ac.classify_axis_status(row))[0]
        regrades[row.variant] = grade
        if row.variant in excluded:
            continue
        if grade != row["mech_consistency_t25"]:
            mismatch.append((row.variant, grade, row["mech_consistency_t25"]))
    check("graded BRCT rows whose regrade differs from stored", len(mismatch), 0)
    check("BRCT partial grades (0.5) present", 
          int((df.loc[brct, "mech_consistency_t25"] == "partial").sum()), 0)
    print("        -> the shipped weighted rubric reproduces every stored grade;")
    print("           the absence of 0.5 is structural (a one-axis variant has no")
    print("           second axis to false-fire), not a different convention.")
    print("           claim is FALSE")

    print("\n[5] REVIEW CLAIM: BRCT lacks uncertainty-gating fields")
    for col in ("ddg_monomer", "ddg_monomer_confident",
                "ddg_monomer_sd", "expected_ddg_monomer"):
        check(f"{col} populated on BRCT", int(df.loc[brct, col].notna().sum()), 12)
    check("BRCT monomer axes entering the denominator", t_brct["monomer"][1], 11)
    print("        -> fields present and evaluated; claim is FALSE")

    print("\n[6] Measured-energy accounting (comparisons, not variants)")
    cal = pd.read_csv(CALIB)
    inb = cal.in_benchmark.astype(str).str.lower().isin(("true", "yes"))
    check("total predicted-measured comparisons", len(cal), 63)
    check("in-benchmark comparisons", int(inb.sum()), 15)
    check("external SKEMPI comparisons", int((~inb).sum()), 48)
    d = cal.measured_kcal.abs() >= 1.0
    check("measured destabilizers >= 1 kcal/mol", int(d.sum()), 44)
    check("of those, in-benchmark", int((d & inb).sum()), 13)
    check("distinct in-benchmark variants measured", int(cal[inb].variant.nunique()), 15)

    print("\n[7] Denominators do NOT all coincide at 57")
    dens = {"MC (variants)": int(graded.sum()),
            "SA (axes, all rows)": sum(v[1] for v in t_all.values()),
            "tier-carrying (variants)": int(df.comavi_tier.notna().sum()),
            "measured (comparisons)": len(cal)}
    print("        " + "  |  ".join(f"{k} = {v}" for k, v in dens.items()))
    check("distinct denominator values", len(set(dens.values())), 4)

    print("\n[8] R1699L/Q sensitivity (ungraded by curated role, not by rubric)")
    check("excluded by curated role", sorted(excluded), ["R1699L", "R1699Q"])
    check("rubric alone would grade them",
          [regrades["R1699L"], regrades["R1699Q"]], ["consistent", "consistent"])
    pts = df.loc[graded, "mech_consistency_t25"].map(GMAP).sum() + 2.0
    check("MC if both included as structurally silent",
          round(pts / 59, 4), 0.7288)
    print("        -> the exclusion is a curation decision applied upstream of the")
    print("           rubric. It must be stated in Methods, not left implicit.")

    if fail:
        print(f"\n{len(fail)} FAILURE(S):")
        for f in fail:
            print("  " + f)
        return 1
    print("\nPASS - every denominator reconciles and both review claims are refuted.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
