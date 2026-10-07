#!/usr/bin/env python3
"""Is there usable signal BELOW the 1.0 kcal/mol screening default?

WHY THIS EXISTS
---------------
A real effect can be small. A variant that weakens an interface by half a
kilocalorie is biologically possible and would fall under every threshold in
the prespecified range, so the paper's calls would be silent on it. The
manuscript says a silent result is not evidence of benignity, but it never
tested the sharper question a reader will ask next: when the predicted energy
is below the screening default, does its VALUE still carry information? If the
sub-threshold energies ranked truly disrupted axes above intact ones, a user
could read them as weak evidence. If they do not, they must be read as no
evidence, which is a different instruction.

This is answerable two ways with what is already committed, and both are
reported because they fail in different places.

AXIS LEVEL -- against the curated ground truth
    Among committed axes whose predicted magnitude falls below the screening
    default, does magnitude separate the destabilizing commitments from the
    neutral ones? Reported as an AUC with a system-cluster bootstrap, because
    the axes come from 14 clustered systems.

MEASURED LEVEL -- against direct measurements
    Among the comparisons with a measured energy, how often does the predicted
    SIGN agree, stratified by predicted magnitude? And how many measured
    effects large enough to matter (1.0 to 2.5 kcal/mol) are predicted below
    the screening default, and so missed?

WHAT THE ANSWER IS USED FOR
---------------------------
The result is a bounded negative, and the manuscript states it as one. The
point estimates lean the right way but the intervals include chance, so the
instruction to a user is that a sub-threshold energy is not weak evidence of
disruption -- it is no evidence either way. That is a stronger and more useful
statement than leaving the question open, and it is the reason the paper does
not offer a sub-threshold reading of its own outputs.

CAVEAT RECORDED IN THE OUTPUT
-----------------------------
The measured set is dominated by two SKEMPI complexes, so the per-band sign
rates rest on few independent systems and no interval is claimed for them.

OUTPUT
------
reference_outputs/COMAVI_subthreshold_evidence.json
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
CALIB = REPO / "reference_outputs" / "COMAVI_delta_calibration_points.csv"
OUT = REPO / "reference_outputs" / "COMAVI_subthreshold_evidence.json"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
SCREENING_DEFAULT = 1.0
REFERENCE_T = 2.5
N_BOOT = 4000
SEED = 11

AXES = {"monomer": ("monomer", "expected_ddg_monomer"),
        "fold_complex": ("fold", "expected_ddg_fold_complex"),
        "binding": ("binding", "expected_ddg_binding")}


def strongest_signed(row, kind, partners):
    """Signed predicted energy of largest magnitude on an axis."""
    if kind == "monomer":
        v = row.get("ddg_monomer")
        return float(v) if pd.notna(v) else None
    vals = [float(row["ddg_%s_%s" % (kind, p)]) for p in partners
            if pd.notna(row.get("ddg_%s_%s" % (kind, p)))]
    return max(vals, key=abs) if vals else None


def axis_table(canon):
    partners = ac.discover_partners(canon)
    graded = canon[canon.mech_consistency_t25.isin(GRADE_MAP)]
    rows = []
    for _, r in graded.iterrows():
        for axis, (kind, col) in AXES.items():
            token = str(r.get(col, "")).strip().lower()
            if token not in ("destab", "stab", "neutral"):
                continue
            v = strongest_signed(r, kind, partners)
            if v is None:
                continue
            rows.append({"system": r.system, "variant": r.variant, "axis": axis,
                         "token": token, "predicted": v, "magnitude": abs(v)})
    return pd.DataFrame(rows)


def auc(df):
    """P(magnitude of a destabilizing axis > that of a neutral axis)."""
    pos = df.loc[df.token == "destab", "magnitude"].to_numpy()
    neg = df.loc[df.token == "neutral", "magnitude"].to_numpy()
    if not len(pos) or not len(neg):
        return None
    return float(stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
                 / (len(pos) * len(neg)))


def build():
    canon = pd.read_csv(CANON, low_memory=False)
    ax = axis_table(canon)
    binary = ax[ax.token.isin(("destab", "neutral"))]
    band = binary[binary.magnitude < SCREENING_DEFAULT]

    rng = np.random.default_rng(SEED)
    systems = sorted(binary.system.unique())
    draws = []
    for _ in range(N_BOOT):
        pick = rng.choice(systems, size=len(systems), replace=True)
        d = pd.concat([binary[binary.system == s] for s in pick], ignore_index=True)
        d = d[d.magnitude < SCREENING_DEFAULT]
        if (d.token == "destab").sum() >= 3 and (d.token == "neutral").sum() >= 3:
            draws.append(auc(d))
    draws = np.array([d for d in draws if d is not None])

    cal = pd.read_csv(CALIB)
    agree = np.sign(cal.foldx_ddg) == np.sign(cal.measured_kcal)
    edges = [0.0, 0.5, SCREENING_DEFAULT, REFERENCE_T, np.inf]
    labels = ["<0.5", "0.5-1.0", "1.0-2.5", ">=2.5"]
    cal["band"] = pd.cut(cal.foldx_ddg.abs(), edges, labels=labels, right=False)
    by_band = [{"predicted_magnitude": lab,
                "n": int((cal.band == lab).sum()),
                "sign_agree": int(agree[cal.band == lab].sum())}
               for lab in labels]

    small = cal[cal.measured_kcal.abs() < SCREENING_DEFAULT]
    mid = cal[(cal.measured_kcal.abs() >= SCREENING_DEFAULT)
              & (cal.measured_kcal.abs() < REFERENCE_T)]

    out = {
        "question": ("When a predicted energy falls below the 1.0 kcal/mol "
                     "screening default, does its value still carry "
                     "information?"),
        "screening_default": SCREENING_DEFAULT,
        "axis_level": {
            "committed_axes_with_a_prediction": int(len(binary)),
            "below_screening_default": int(len(band)),
            "below_and_committed_destabilizing": int((band.token == "destab").sum()),
            "below_and_committed_neutral": int((band.token == "neutral").sum()),
            "auc_full_range": round(auc(binary), 3),
            "auc_within_band": round(auc(band), 3),
            "auc_within_band_ci95": [round(float(np.percentile(draws, 2.5)), 3),
                                     round(float(np.percentile(draws, 97.5)), 3)],
            "p_auc_at_or_below_chance": round(float((draws <= 0.5).mean()), 3),
            "median_magnitude_destab": round(float(band.loc[band.token == "destab",
                                                            "magnitude"].median()), 3),
            "median_magnitude_neutral": round(float(band.loc[band.token == "neutral",
                                                             "magnitude"].median()), 3),
            "n_bootstrap_draws": int(len(draws)),
        },
        "measured_level": {
            "n_comparisons": int(len(cal)),
            "sign_agreement_by_predicted_band": by_band,
            "measured_below_default": {
                "n": int(len(small)),
                "also_predicted_below_default": int((small.foldx_ddg.abs()
                                                     < SCREENING_DEFAULT).sum()),
                "sign_agree": int((np.sign(small.foldx_ddg)
                                   == np.sign(small.measured_kcal)).sum()),
            },
            "measured_between_default_and_reference": {
                "n": int(len(mid)),
                "predicted_below_default": int((mid.foldx_ddg.abs()
                                                < SCREENING_DEFAULT).sum()),
                "note": ("Measured effects large enough to matter that the "
                         "screening default does not reach."),
            },
            "caveat": ("48 of the 63 comparisons come from two SKEMPI "
                       "complexes, so the per-band rates rest on few "
                       "independent systems; no interval is claimed for them."),
        },
        "conclusion": ("The point estimates lean the right way but the "
                       "interval on the within-band AUC includes chance, so a "
                       "sub-threshold energy is not weak evidence of "
                       "disruption; it is no evidence either way."),
    }

    # The manuscript states this as a bounded negative. If a future change
    # makes the band informative, the claim must be rewritten rather than
    # silently kept, so assert the shape the sentence depends on.
    lo, hi = out["axis_level"]["auc_within_band_ci95"]
    assert lo <= 0.5 <= hi, ("the within-band AUC interval no longer includes "
                             "chance; the Limitations sentence must be rewritten")
    assert out["axis_level"]["below_and_committed_destabilizing"] > 0, \
        "no committed destabilizing axis falls below the default; claim is moot"
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
            print("FAIL: sub-threshold record is stale; differs at: %s"
                  % [k for k in fresh if prev.get(k) != fresh.get(k)])
            return 1
        print("PASS: sub-threshold evidence reproduces "
              "(within-band AUC %.3f, 95%% CI %s)"
              % (fresh["axis_level"]["auc_within_band"],
                 fresh["axis_level"]["auc_within_band_ci95"]))
        return 0
    OUT.write_text(json.dumps(fresh, indent=1) + "\n")
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
