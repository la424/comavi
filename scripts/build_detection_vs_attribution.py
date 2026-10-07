#!/usr/bin/env python3
"""Does COMAVI notice a lesion, and does it name the right one?

WHY THIS EXISTS
---------------
The headline mechanism-pattern score answers neither question on its own. It is
a partial-credit rubric over a cohort that is half curated-silent variants, so
it sits between two quantities a reader actually wants:

  detection    on a variant with a curated structural mechanism, does ANY of
               the three axes fire at the reference threshold?
  attribution  does the axis the mechanism is attributed to fire?

Separating them is what shows where the pipeline is strong and where it is not,
and it localises the weakest arm to a specific failure rather than a low number.

DEFINITIONS
-----------
Firing is |ddG| >= the reference threshold on an axis, using the monomer value
and the strongest-magnitude partner for the complex-context fold and binding
axes -- the same reduction the concordance layer uses.

The expected axis follows the curated mechanism class:
  fold_mechanism        monomer fold
  ppi_destab_mechanism  binding
  mixed_structural      at least one of {complex fold, binding} AND at least
                        one of {monomer, complex fold} -- a mixed expectation
                        is not satisfied by a single axis

Correct rejection is the complement on curated structurally silent variants:
no axis fires.

Usage
-----
    PYTHONPATH=. python scripts/build_detection_vs_attribution.py
    PYTHONPATH=. python scripts/build_detection_vs_attribution.py --check
"""
import argparse
import json
import pathlib

import pandas as pd

# The axis-status vocabulary is owned by the concordance helper. Importing it
# rather than re-deriving "which axis does the literature commit" keeps this
# module from drifting away from the grading path, which is how a hand-mirrored
# copy went stale earlier in this project.
from apply_concordance_v5 import classify_axis_status

REPO = pathlib.Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_detection_vs_attribution.json"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
REFERENCE_T = 2.5
# v8.0 adds "ppi_stab_mechanism". Until SMAD4 I500V there was no stabilizing
# commitment in the graded set, so the tuple was complete by accident rather
# than by design, and a stabilizing expectation fell through both buckets --
# which the population identity below caught. A strengthened interaction is a
# curated structural lesion like any other, and detection asks only whether an
# axis fired, which is tested on |ddG|.
STRUCTURAL = ("fold_mechanism", "ppi_destab_mechanism", "mixed_structural",
              "ppi_stab_mechanism")
SILENT = "structurally_silent"


def axis_columns(canon, kind):
    # "_indistinguishable" is listed explicitly. The original exclusion was
    # "_disting", which does NOT match "ddg_binding_<partner>_indistinguishable"
    # -- after that underscore comes "i", not "d". So all 25 per-partner
    # indistinguishability FLAGS were being treated as binding energies, and
    # because strongest() takes the maximum by absolute value, a flag of 1.0 won
    # whenever no real partner energy exceeded 1.0 kcal/mol. That fired a
    # spurious binding call at the 1.0 threshold and nowhere else, which is why
    # correct rejection at 1.0 read 5/32 while the rubric arm -- computed by
    # apply_concordance_v5, which excludes the flags correctly -- read 14/32 on
    # the same property. Two implementations of one quantity disagreeing by nine
    # variants was the visible symptom; this was the cause.
    #
    # This module is the only one that builds a ddG column list from scratch.
    # Every other caller starts from apply_concordance_v5.discover_partners(),
    # which already excludes the flags, so the defect did not spread.
    bad = ("_sd", "_ci95", "_disting", "_indistinguishable",
           "_runs", "_vote", "_confident")
    pre = "ddg_%s_" % kind
    cols = [c for c in canon.columns
            if c.startswith(pre) and not any(b in c for b in bad)]
    leaked = [c for c in cols if canon[c].dropna().isin([0, 1, True, False]).all()
              and canon[c].notna().sum() > 2 * len(canon) / 3]
    assert not leaked, ("these look like flags, not energies: %s" % leaked[:4])
    return cols


def strongest(row, cols, prefix):
    vals = [row[c] for c in cols if pd.notna(row.get(c))]
    return max(vals, key=abs) if vals else None


def per_axis_thresholds(threshold):
    """Expand a threshold spec into (monomer, fold, binding).

    Four of the five committed thresholds are scalar; the substitution-adjusted
    one is per-axis ({'monomer': 2.9, 'fold': 2.9, 'binding': 3.5}). Accepting
    both here is what lets build() serve the whole committed set, so the
    threshold sweep is this same computation evaluated at five points rather
    than a second implementation that agrees at one.
    """
    if isinstance(threshold, dict):
        return (float(threshold["monomer"]), float(threshold["fold"]),
                float(threshold["binding"]))
    t = float(threshold)
    return t, t, t


def build(canon, threshold=REFERENCE_T):
    graded = canon[canon.mech_consistency_t25.isin(GRADE_MAP)].copy()
    graded["score"] = graded.mech_consistency_t25.map(GRADE_MAP)
    fold_cols = axis_columns(canon, "fold")
    bind_cols = axis_columns(canon, "binding")
    graded["fold_max"] = [strongest(r, fold_cols, "ddg_fold_") for _, r in graded.iterrows()]
    graded["bind_max"] = [strongest(r, bind_cols, "ddg_binding_") for _, r in graded.iterrows()]

    t_mono, t_fold, t_bind = per_axis_thresholds(threshold)

    def fires_at(t):
        def f(v):
            return bool(pd.notna(v) and abs(float(v)) >= t)
        return f

    graded["mono_fires"] = graded.ddg_monomer.map(fires_at(t_mono))
    graded["fold_fires"] = graded.fold_max.map(fires_at(t_fold))
    graded["bind_fires"] = graded.bind_max.map(fires_at(t_bind))
    graded["any_fires"] = graded[["mono_fires", "fold_fires", "bind_fires"]].any(axis=1)

    # v8.6: fold attribution is now decided by WHICH fold axis the literature
    # commits, not by assuming the isolated-subunit one.
    #
    # derive_expected_mech_class puts a variant in fold_mechanism when EITHER
    # fold axis is committed positive -- isolated-subunit or complex-context --
    # but this rule asked only `mono_fires`. A lesion committed on the
    # complex-context fold axis therefore could not be credited no matter what
    # the pipeline predicted, so the middle axis of a three-axis method had no
    # route to a correct attribution. BRCA1 C61G is the case: committed on both
    # fold axes, isolated-subunit prediction 2.33 (below the 2.5 threshold) and
    # complex-context 6.78 (well above). It was scored as a wrong attribution
    # while naming exactly the axis the literature names.
    #
    # The fix credits a fold mechanism when a COMMITTED fold axis fires. It is
    # not a loosening: an uncommitted axis still cannot earn credit, and
    # MLH1 R755W -- committed on both fold axes, predicted 0.97 and 0.94 with
    # binding at 9.66 -- remains a miss, because the pipeline names binding
    # where the literature names fold. Attribution moves 14 -> 15 of 29.
    axis_state = {(r.system, r.variant): classify_axis_status(r)
                  for _, r in graded.iterrows()}

    def committed_positive(r, axis):
        st = axis_state.get((r.system, r.variant)) or {}
        return st.get(axis) == "positive"

    def right_axis(r):
        if r.expected_mech_class == "fold_mechanism":
            hits = []
            if committed_positive(r, "fold_monomer"):
                hits.append(bool(r.mono_fires))
            if committed_positive(r, "fold_complex"):
                hits.append(bool(r.fold_fires))
            # A fold_mechanism always has at least one committed fold axis;
            # assert rather than fall back, so a vocabulary change upstream
            # surfaces here instead of silently scoring every fold variant
            # as a miss.
            assert hits, ("fold_mechanism with no committed fold axis: %s %s"
                          % (r.system, r.variant))
            return any(hits)
        if r.expected_mech_class == "ppi_destab_mechanism":
            return bool(r.bind_fires)
        if r.expected_mech_class == "mixed_structural":
            return bool(r.fold_fires or r.bind_fires) and bool(r.mono_fires or r.fold_fires)
        return None

    graded["right_axis"] = graded.apply(right_axis, axis=1)
    struct = graded[graded.expected_mech_class.isin(STRUCTURAL)]
    silent = graded[graded.expected_mech_class == SILENT]

    # What scoring the assembled complex adds, and what it costs. A monomer-only
    # analysis sees only the isolated-subunit axis, so detection and correct
    # rejection are recounted on that axis alone, on the same variants and at the
    # same threshold. "In a complex" means at least one partner axis was computed;
    # for the single-chain BRCT system the two counts coincide by construction.
    in_complex = graded.fold_max.notna() | graded.bind_max.notna()
    struct_cx = struct[in_complex.loc[struct.index]]
    only_assembled = struct[struct.any_fires & ~struct.mono_fires]

    per_class = {}
    for cls, sub in graded.groupby("expected_mech_class"):
        row = {"n": int(len(sub)),
               "mean_grade": round(float(sub.score.mean()), 4),
               "consistent": int((sub.mech_consistency_t25 == "consistent").sum()),
               "partial": int((sub.mech_consistency_t25 == "partial").sum()),
               "inconsistent": int((sub.mech_consistency_t25 == "inconsistent").sum())}
        if cls in STRUCTURAL:
            row["detection"] = int(sub.any_fires.sum())
            row["attribution"] = int(sub.right_axis.sum())
            # v8.6: separate accounting for the axis each commitment names, so
            # the complex-context fold axis is visible rather than folded into
            # an isolated-subunit total. `n` is how many variants in the class
            # commit that axis; `fires` is how many of those have it fire.
            by_axis = {}
            for label, axis, col in (("isolated_subunit", "fold_monomer", "mono_fires"),
                                     ("complex_fold", "fold_complex", "fold_fires"),
                                     ("binding", "binding", "bind_fires")):
                cm = sub.apply(lambda r: committed_positive(r, axis), axis=1)
                if int(cm.sum()):
                    by_axis[label] = {"n": int(cm.sum()),
                                      "fires": int(sub.loc[cm, col].sum())}
            row["committed_axis_accounting"] = by_axis
        elif cls == SILENT:
            row["correct_rejection"] = int((~sub.any_fires).sum())
        per_class[cls] = row

    ppi = graded[graded.expected_mech_class == "ppi_destab_mechanism"]
    return {
        "threshold": threshold,
        "population": {"graded": int(len(graded)),
                       "structural": int(len(struct)),
                       "silent": int(len(silent))},
        "headline_score": {"total": float(graded.score.sum()),
                           "n": int(len(graded)),
                           "mean": round(float(graded.score.mean()), 4)},
        "detection": {"k": int(struct.any_fires.sum()), "n": int(len(struct)),
                      "rate": round(float(struct.any_fires.mean()), 3)},
        "attribution": {"k": int(struct.right_axis.sum()), "n": int(len(struct)),
                        "rate": round(float(struct.right_axis.mean()), 3)},
        "correct_rejection": {"k": int((~silent.any_fires).sum()), "n": int(len(silent)),
                              "rate": round(float((~silent.any_fires).mean()), 3)},
        "isolated_subunit_only": {
            "detection": {"k": int(struct.mono_fires.sum()), "n": int(len(struct)),
                          "rate": round(float(struct.mono_fires.mean()), 3)},
            "correct_rejection": {"k": int((~silent.mono_fires).sum()), "n": int(len(silent)),
                                  "rate": round(float((~silent.mono_fires).mean()), 3)},
            "detected_only_with_assembled_axes": {
                "k": int(len(only_assembled)),
                "variants": ["%s %s" % (str(r.gene).upper(), r.variant)
                             for _, r in only_assembled.iterrows()]},
            "structural_in_complexes": {
                "n": int(len(struct_cx)),
                "detection_all_axes": int(struct_cx.any_fires.sum()),
                "detection_isolated_subunit": int(struct_cx.mono_fires.sum())},
        },
        "per_class": per_class,
        "interface_class_diagnosis": {
            "n": int(len(ppi)),
            "binding_axis_fires": int(ppi.bind_fires.sum()),
            "complex_fold_axis_fires": int(ppi.fold_fires.sum()),
            "note": ("For variants curated as pure interface destabilizers the "
                     "complex-context fold axis fires more often than the binding "
                     "axis, so the lesion is detected but attributed to the wrong "
                     "axis. This is the separation the pipeline exists to make and "
                     "it is where the benchmark is weakest."),
            "large_monomer_effect": [
                {"gene": str(r.gene).upper(), "variant": r.variant,
                 "ddg_monomer": round(float(r.ddg_monomer), 2)}
                for _, r in ppi.iterrows()
                if pd.notna(r.ddg_monomer) and abs(float(r.ddg_monomer)) >= 5.0],
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    result = build(canon)

    # Structural checks, not frozen literals. The v7.9 ground-truth correction
    # moved the headline from 39.5 to 41.0 and three literals like the one that
    # used to sit here all failed at once; bumping them would have removed the
    # guard. Recomputing the expectation from the canonical keeps it live.
    _W = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
    _g = canon[canon.mech_consistency_t25.isin(_W)]
    p = result["population"]
    assert p["graded"] == len(_g), (p["graded"], len(_g))
    assert p["structural"] + p["silent"] == p["graded"]
    assert abs(result["headline_score"]["total"]
               - _g.mech_consistency_t25.map(_W).sum()) < 1e-9
    # The weakest-arm claim in the manuscript rests on this ordering. Restrict
    # it to classes with at least MIN_ARM_N members: v8.0 introduced
    # ppi_stab_mechanism with a single variant, and a one-member class scoring
    # zero would take the title and make the manuscript assert a ranking that
    # n = 1 cannot support. The floor is named rather than implicit so that a
    # class crossing it changes the claim deliberately.
    MIN_ARM_N = 3
    means = {k: v["mean_grade"] for k, v in result["per_class"].items()
             if v["n"] >= MIN_ARM_N}
    assert min(means, key=means.get) == "ppi_destab_mechanism", means
    assert all(v["n"] < MIN_ARM_N or k in STRUCTURAL + (SILENT,)
               for k, v in result["per_class"].items()), "unbucketed class"
    # The isolated-subunit recount is a subset of the full profile: everything it
    # detects the full profile detects, and the remainder is exactly the set
    # detected only through the assembled-complex and binding axes.
    iso = result["isolated_subunit_only"]
    assert (iso["detection"]["k"] + iso["detected_only_with_assembled_axes"]["k"]
            == result["detection"]["k"]), iso
    assert iso["correct_rejection"]["k"] >= result["correct_rejection"]["k"], iso

    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: detection/attribution matches a fresh build")
        return 0

    print("graded %d (structural %d, silent %d) at threshold %.1f kcal/mol"
          % (p["graded"], p["structural"], p["silent"], result["threshold"]))
    for k in ("detection", "attribution", "correct_rejection"):
        d = result[k]
        print("  %-18s %2d/%-2d = %.3f" % (k, d["k"], d["n"], d["rate"]))
    iso = result["isolated_subunit_only"]
    print("  isolated subunit only: detection %d/%d, correct rejection %d/%d; "
          "%d detected only through the assembled axes"
          % (iso["detection"]["k"], iso["detection"]["n"], iso["correct_rejection"]["k"],
             iso["correct_rejection"]["n"], iso["detected_only_with_assembled_axes"]["k"]))
    print("\nper expected mechanism class:")
    for cls, r in sorted(result["per_class"].items(),
                         key=lambda kv: -kv[1]["mean_grade"]):
        extra = ""
        if "detection" in r:
            extra = "detect %2d/%-2d  attribute %2d/%-2d" % (
                r["detection"], r["n"], r["attribution"], r["n"])
        elif "correct_rejection" in r:
            extra = "reject %2d/%-2d" % (r["correct_rejection"], r["n"])
        print("  %-22s n=%2d  grade %.3f   %s" % (cls, r["n"], r["mean_grade"], extra))
    d = result["interface_class_diagnosis"]
    print("\ninterface class (n=%d): binding axis fires %d, complex-fold axis fires %d"
          % (d["n"], d["binding_axis_fires"], d["complex_fold_axis_fires"]))
    if d["large_monomer_effect"]:
        print("  large monomer effects among them: %s"
              % ", ".join("%s %s (%.1f)" % (x["gene"], x["variant"], x["ddg_monomer"])
                          for x in d["large_monomer_effect"]))

    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
