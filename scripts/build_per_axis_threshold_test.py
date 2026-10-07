#!/usr/bin/env python3
"""Do per-axis thresholds recover the interface attribution shortfall?

WHY THIS EXISTS
---------------
The Discussion called per-axis thresholds "the most promising single change to
the calling layer" and then said whether they recover part of the interface
shortfall "has not been tested". A reviewer asked why not. It is testable with
what is already in the repository, and leaving it untested understates a
limitation rather than overstating a result.

WHAT IS TESTED
--------------
Binding agreement peaks at the permissive end of the prespecified range while
the two stability axes peak at the stringent end, so the obvious per-axis
scheme holds the stability axes at the 2.5 kcal/mol reference and lowers the
binding threshold. This sweeps the binding threshold across the prespecified
range and reports detection, attribution, interface attribution and correct
rejection at each point.

THE ANSWER, AND WHY IT IS A LIMITATION RATHER THAN A FIX
--------------------------------------------------------
Lowering the binding threshold does recover part of the shortfall -- but only
part, and the reason the rest is unreachable is the point. Half the interface
mechanisms carry a predicted binding energy below the most permissive
threshold in the prespecified range, so NO placement of a binding threshold
reaches them. For those variants the calling layer is not what fails; the
predicted energy is near zero where the literature reports a disrupted
interface. That localises the shortfall to the energy prediction, which is a
sharper statement than "per-axis thresholds might help".

IN-SAMPLE
---------
The recovered operating point is selected on the same 56 variants it is
evaluated on, so the gain is an upper bound, not an out-of-sample estimate.
The unreachable fraction is NOT subject to that caveat: it is a property of
the predicted energies, independent of which threshold is chosen.

OUTPUT
------
reference_outputs/COMAVI_per_axis_threshold_test.json
"""
import argparse
import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_per_axis_threshold_test.json"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
REFERENCE_T = 2.5
# The prespecified range, fixed before any variant was scored.
PRESPECIFIED = [1.0, 1.5, 2.0, 2.5]
INTERFACE_CLASSES = ("ppi_destab_mechanism", "ppi_stab_mechanism")
SILENT_CLASS = "structurally_silent"


def axis_magnitude(row, kind, partners):
    """Largest absolute predicted energy on an axis, 0.0 when the axis is absent."""
    if kind == "monomer":
        v = row.get("ddg_monomer")
        return abs(float(v)) if pd.notna(v) else 0.0
    vals = [abs(float(row["ddg_%s_%s" % (kind, p)])) for p in partners
            if pd.notna(row.get("ddg_%s_%s" % (kind, p)))]
    return max(vals) if vals else 0.0


def build_frame(canon):
    partners = ac.discover_partners(canon)
    graded = canon[canon.mech_consistency_t25.isin(GRADE_MAP)]
    rows = []
    for _, r in graded.iterrows():
        state = ac.classify_axis_status(r)
        rows.append({
            "system": r.system, "variant": r.variant,
            "mech_class": r.expected_mech_class,
            "monomer": axis_magnitude(r, "monomer", partners),
            "complex_fold": axis_magnitude(r, "fold", partners),
            "binding": axis_magnitude(r, "binding", partners),
            "commits_monomer": state["fold_monomer"] == "positive",
            "commits_complex_fold": state["fold_complex"] == "positive",
            "commits_binding": state["binding"] == "positive",
        })
    df = pd.DataFrame(rows)
    df["is_silent"] = df.mech_class == SILENT_CLASS
    df["is_interface"] = df.mech_class.isin(INTERFACE_CLASSES)
    return df


def evaluate(df, t_monomer, t_complex_fold, t_binding):
    """Detection, attribution, interface attribution and correct rejection.

    Attribution follows the same rule the shipped detection generator uses: a
    fold mechanism is credited when a COMMITTED fold axis fires, an interface
    mechanism when the binding axis fires, and a mixed mechanism when an
    assembled axis and a fold axis both fire.
    """
    fires_m = df.monomer >= t_monomer
    fires_f = df.complex_fold >= t_complex_fold
    fires_b = df.binding >= t_binding
    fires_any = fires_m | fires_f | fires_b
    structural = ~df.is_silent

    attributed = (
        ((df.mech_class == "fold_mechanism")
         & ((df.commits_monomer & fires_m) | (df.commits_complex_fold & fires_f)))
        | (df.is_interface & fires_b)
        | ((df.mech_class == "mixed_structural") & (fires_f | fires_b) & (fires_m | fires_f))
    )
    return {
        "t_monomer": t_monomer, "t_complex_fold": t_complex_fold, "t_binding": t_binding,
        "detection_k": int((fires_any & structural).sum()),
        "detection_n": int(structural.sum()),
        "attribution_k": int(attributed.sum()),
        "attribution_n": int(structural.sum()),
        "interface_attribution_k": int((df.is_interface & fires_b).sum()),
        "interface_attribution_n": int(df.is_interface.sum()),
        "correct_rejection_k": int((~fires_any & df.is_silent).sum()),
        "correct_rejection_n": int(df.is_silent.sum()),
    }


def build():
    canon = pd.read_csv(CANON, low_memory=False)
    df = build_frame(canon)

    reference = evaluate(df, REFERENCE_T, REFERENCE_T, REFERENCE_T)
    sweep = [evaluate(df, REFERENCE_T, REFERENCE_T, tb) for tb in PRESPECIFIED]

    # Best binding threshold by interface attribution; ties go to the most
    # stringent, so the reported gain is never inflated by a looser tie.
    best = max(sweep, key=lambda r: (r["interface_attribution_k"], r["t_binding"]))

    iface = df[df.is_interface]
    floor = min(PRESPECIFIED)
    unreachable = iface[iface.binding < floor]
    # Of the unreachable ones, how many does the complex-fold axis still catch?
    caught_elsewhere = int((unreachable.complex_fold >= REFERENCE_T).sum())

    out = {
        "question": ("Do per-axis thresholds recover the interface attribution "
                     "shortfall left at the uniform reference threshold?"),
        "reference_uniform": reference,
        "binding_threshold_sweep": sweep,
        "best_by_interface_attribution": best,
        "gain_over_reference": {
            "attribution": best["attribution_k"] - reference["attribution_k"],
            "interface_attribution": (best["interface_attribution_k"]
                                      - reference["interface_attribution_k"]),
            "correct_rejection": (best["correct_rejection_k"]
                                  - reference["correct_rejection_k"]),
            "in_sample": True,
            "note": ("The operating point is selected on the variants it is "
                     "evaluated on, so the gain is an upper bound."),
        },
        "unreachable_by_any_threshold": {
            "k": int(len(unreachable)),
            "n": int(len(iface)),
            "floor_kcal_mol": floor,
            "median_predicted_binding": round(float(unreachable.binding.median()), 3),
            "variants": sorted("%s %s" % (r.system, r.variant)
                               for _, r in unreachable.iterrows()),
            "detected_on_complex_fold_at_reference": caught_elsewhere,
            "note": ("Predicted binding energy below the most permissive "
                     "prespecified threshold, so no placement of a binding "
                     "threshold reaches them. Not an in-sample quantity."),
        },
    }

    # The conclusion the manuscript draws rests on these two being true
    # together: thresholds help, and they cannot close the gap.
    assert out["gain_over_reference"]["interface_attribution"] > 0, \
        "no per-axis gain; the manuscript sentence would be wrong"
    assert len(unreachable) * 2 >= len(iface), \
        "unreachable fraction is no longer at least half; rewrite the claim"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="recompute and fail if the committed record differs")
    args = ap.parse_args()
    fresh = json.loads(json.dumps(build()))
    if args.check:
        if not OUT.exists():
            print("FAIL: %s is missing" % OUT.name)
            return 1
        prev = json.loads(OUT.read_text())
        if prev != fresh:
            moved = [k for k in fresh if prev.get(k) != fresh.get(k)]
            print("FAIL: per-axis threshold record is stale; differs at: %s" % moved)
            return 1
        print("PASS: per-axis threshold test reproduces "
              "(interface attribution %d/%d at the reference, %d unreachable)"
              % (fresh["reference_uniform"]["interface_attribution_k"],
                 fresh["reference_uniform"]["interface_attribution_n"],
                 fresh["unreachable_by_any_threshold"]["k"]))
        return 0
    OUT.write_text(json.dumps(fresh, indent=1) + "\n")
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
