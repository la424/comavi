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
The obvious per-axis scheme holds the two stability axes at the 2.5 kcal/mol
reference and lowers the binding threshold. This sweeps the binding threshold
across the prespecified range and reports detection, attribution, interface
attribution and correct rejection at each point, and lists the interface
mechanisms still missed at the reference with their predicted energies.

THE ANSWER (v11)
----------------
Before v11 the binding threshold appeared to matter: lowering it to 1.5 kcal/mol
raised interface attribution from 4 to 6 of 12. Both gains were BRCA1 V11G and
MLH1 L749P, whose evidence measures an amount of complex. Under evidence
resolution the assembled-complex stability axis already matches them at the
reference, so lowering the binding threshold recovers nothing and costs correct
rejections. Every interface mechanism still missed at the reference carries a
predicted binding energy below the most permissive prespecified threshold, so
NO placement of a binding threshold reaches it. The calling layer is not what
fails there: the predicted energy is near zero where the literature reports a
changed interaction, which locates the shortfall in the energy prediction.

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
import importlib.util
import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import apply_concordance_v5 as ac  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "build_detection_vs_attribution", REPO / "scripts" / "build_detection_vs_attribution.py")
bda = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bda)

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_per_axis_threshold_test.json"
DETECTION_RECORD = REPO / "reference_outputs" / "COMAVI_detection_vs_attribution.json"

# v11 NOTE. Until v11 this file restated the attribution rule instead of
# calling the detection generator, so when evidence resolution changed that rule
# (an interface expectation measured only as an amount of complex is matched by
# either assembled axis), this record kept reporting attribution 15/29 at the
# reference while the headline record reported 18/29. Rates now come from
# build_detection_vs_attribution.build() at a per-axis threshold, so there is one
# implementation, and the reference point is checked against the headline record.

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


def evaluate(canon, t_monomer, t_complex_fold, t_binding):
    """Detection, attribution, interface attribution and correct rejection at a
    per-axis threshold, from the shipped detection generator (one rule for
    attribution, evidence resolution included). Returns (rates, frame)."""
    r, frame = bda.build(canon, {"monomer": t_monomer, "fold": t_complex_fold,
                                 "binding": t_binding}, return_frame=True)
    iface = frame.expected_mech_class.isin(INTERFACE_CLASSES)
    right = frame.right_axis.map(lambda v: v is True)
    rates = {
        "t_monomer": t_monomer, "t_complex_fold": t_complex_fold, "t_binding": t_binding,
        "detection_k": int(r["detection"]["k"]), "detection_n": int(r["detection"]["n"]),
        "attribution_k": int(r["attribution"]["k"]), "attribution_n": int(r["attribution"]["n"]),
        "interface_attribution_k": int((iface & right).sum()),
        "interface_attribution_n": int(iface.sum()),
        "correct_rejection_k": int(r["correct_rejection"]["k"]),
        "correct_rejection_n": int(r["correct_rejection"]["n"]),
    }
    return rates, frame


def build():
    canon = pd.read_csv(CANON, low_memory=False)
    df = build_frame(canon)

    reference, ref_frame = evaluate(canon, REFERENCE_T, REFERENCE_T, REFERENCE_T)
    headline = json.loads(DETECTION_RECORD.read_text())
    assert (reference["detection_k"], reference["attribution_k"], reference["correct_rejection_k"]) == (
        headline["detection"]["k"], headline["attribution"]["k"], headline["correct_rejection"]["k"]), \
        "the reference point must reproduce COMAVI_detection_vs_attribution.json"
    sweep = [evaluate(canon, REFERENCE_T, REFERENCE_T, tb)[0] for tb in PRESPECIFIED]

    # Best binding threshold by interface attribution; ties go to the most
    # stringent, so the reported gain is never inflated by a looser tie.
    best = max(sweep, key=lambda r: (r["interface_attribution_k"], r["t_binding"]))

    # Interface mechanisms not attributed at the reference (undetected, or
    # detected on another axis), and whether any binding threshold in the
    # prespecified range could reach them.
    ref_right = ref_frame.set_index(["system", "variant"]).right_axis.map(lambda v: v is True)
    iface = df[df.is_interface].copy()
    iface["attributed_at_reference"] = [bool(ref_right.loc[(s, v)])
                                        for s, v in zip(iface.system, iface.variant)]
    floor = min(PRESPECIFIED)
    missed = iface[~iface.attributed_at_reference]
    unreachable = missed[missed.binding < floor]
    # Interface mechanisms whose predicted binding is below the floor but that
    # are attributed at the reference anyway, through assembled-complex stability
    # (evidence measured as an amount of complex).
    caught_elsewhere = int((iface.attributed_at_reference & (iface.binding < floor)).sum())

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
        "interface_not_attributed_at_reference": {
            "k": int(len(missed)),
            "n": int(len(iface)),
            "variants": [{"system": r.system, "variant": r.variant,
                          "predicted_binding": round(float(r.binding), 3),
                          "predicted_complex_fold": round(float(r.complex_fold), 3)}
                         for _, r in missed.sort_values(["system", "variant"]).iterrows()],
        },
        "unreachable_by_any_threshold": {
            "k": int(len(unreachable)),
            "of_missed": int(len(missed)),
            "n": int(len(iface)),
            "floor_kcal_mol": floor,
            "median_predicted_binding": (round(float(unreachable.binding.median()), 3)
                                         if len(unreachable) else None),
            "variants": sorted("%s %s" % (r.system, r.variant)
                               for _, r in unreachable.iterrows()),
            "below_floor_but_attributed_via_assembled_stability": caught_elsewhere,
            "note": ("Interface mechanisms not attributed at the reference whose predicted "
                     "binding energy is below the most permissive prespecified "
                     "threshold, so no placement of a binding threshold reaches "
                     "them. Not an in-sample quantity."),
        },
    }
    # Internal consistency only: what the manuscript says about these values is
    # asserted where the sentence is written, from this record.
    assert len(unreachable) <= len(missed) <= len(iface)
    assert out["gain_over_reference"]["interface_attribution"] >= 0
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
