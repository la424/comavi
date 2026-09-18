#!/usr/bin/env python3
"""Directional asymmetry between firing and silent structural outputs.

WHY THIS EXISTS
---------------
S6 Text reported this block entirely from hand-typed numbers: no file in the
repository contained "25/61", "36/61", "3.1 x 10^-4" or "11/14 systems". On
audit, four of its five figures did not reproduce against the released
canonical, while the one that carried the inferential claim -- that every
clinically labelled variant with a firing output was pathogenic -- did. The
stated numbers were internally self-consistent (22 + 20 + 14 = 56), which is
why nothing flagged them: they described a labelled population of 56 that the
canonical does not contain.

WHAT IT COMPUTES
----------------
Firing uses the shipped reduction, apply_concordance_v5.compute_max_abs_ddg,
which counts only confidence-passing axes and drops binding values that are
statistically indistinguishable from zero. Recomputing it here rather than
thresholding raw ddG columns matters: a naive |ddG| >= threshold sweep over all
partner columns gives a different count because it ignores those gates.

Clinical ground truth is the curated `phenotype` annotation, with
gain-of-function counted as pathogenic. This is clinical, NOT the structural
ground truth used for COMAVI's mechanism metrics, so it was unaffected by the
v7.6 ledger correction.

The claim the block exists to support is the DIRECTION of the asymmetry, not a
predictive value: the benchmark was selected for mechanism evidence, so these
fractions are not expected clinical performance.

Usage
-----
    PYTHONPATH=. python scripts/build_directional_enrichment.py
    PYTHONPATH=. python scripts/build_directional_enrichment.py --check
"""
import argparse
import json
import pathlib
import sys

import pandas as pd
from scipy.stats import fisher_exact

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "reference_outputs" / "COMAVI_directional_enrichment.json"
REFERENCE_T = 2.5


def real_partners(canon):
    noise = ("_ci95", "_distinguishable", "vote_strict", "vote_relaxed")
    return sorted({p for p in ac.discover_partners(canon)
                   if p and not any(n in p for n in noise)})


def build(canon, threshold=REFERENCE_T):
    partners = real_partners(canon)
    canon = canon.copy()
    canon["max_abs_ddg"] = canon.apply(
        lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)
    canon["fires"] = canon.max_abs_ddg >= threshold

    lab = canon[canon.phenotype.notna()].copy()
    lab["pathogenic"] = lab.phenotype.astype(str).str.startswith("pathogenic")

    fp = int((lab.fires & lab.pathogenic).sum())
    fb = int((lab.fires & ~lab.pathogenic).sum())
    sp = int((~lab.fires & lab.pathogenic).sum())
    sb = int((~lab.fires & ~lab.pathogenic).sum())
    odds, p = fisher_exact([[fp, fb], [sp, sb]])

    silent_path = lab[~lab.fires & lab.pathogenic]
    return {
        "threshold": threshold,
        "firing_definition": ("apply_concordance_v5.compute_max_abs_ddg >= threshold; "
                              "confidence-passing axes only, indistinguishable binding "
                              "values excluded"),
        "clinical_ground_truth": "curated phenotype; gain-of-function counted as pathogenic",
        "cohort": {"n": int(len(canon)),
                   "firing": int(canon.fires.sum()),
                   "silent": int((~canon.fires).sum())},
        "labelled": {"n": int(len(lab)),
                     "pathogenic": int(lab.pathogenic.sum()),
                     "benign": int((~lab.pathogenic).sum())},
        "table_2x2": {"firing_pathogenic": fp, "firing_benign": fb,
                      "silent_pathogenic": sp, "silent_benign": sb},
        "fisher": {"odds_ratio": (None if odds == float("inf") else round(float(odds), 3)),
                   "odds_ratio_is_infinite": odds == float("inf"),
                   "p_value": float("%.3g" % p)},
        "pathogenic_but_silent": {
            "n": int(len(silent_path)),
            "systems_spanned": int(silent_path.system.nunique()),
            "systems_total": int(canon.system.nunique()),
            "variants": ["%s %s" % (str(r.gene).upper(), r.variant)
                         for _, r in silent_path.iterrows()],
        },
        "interpretation": ("Every clinically labelled variant with a firing structural "
                           "output was pathogenic, so a firing output was never observed "
                           "on a benign variant here; silence carried both classes. The "
                           "asymmetry is directional evidence only -- the cohort was "
                           "selected for mechanism evidence, so these are not expected "
                           "clinical predictive values."),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    result = build(canon)

    t = result["table_2x2"]
    assert t["firing_pathogenic"] + t["firing_benign"] \
        + t["silent_pathogenic"] + t["silent_benign"] == result["labelled"]["n"]
    assert result["cohort"]["firing"] + result["cohort"]["silent"] == result["cohort"]["n"]
    # the inferential claim the block carries
    assert t["firing_benign"] == 0, "a benign variant now has a firing output"

    if args.check:
        if not OUT.exists():
            print("FAIL: %s missing" % OUT.name)
            return 1
        if json.loads(OUT.read_text()) != result:
            print("FAIL: %s differs from a fresh build" % OUT.name)
            return 1
        print("PASS: directional enrichment matches a fresh build")
        return 0

    c, l = result["cohort"], result["labelled"]
    print("cohort: %d firing / %d silent of %d at %.1f kcal/mol"
          % (c["firing"], c["silent"], c["n"], result["threshold"]))
    print("clinically labelled: %d (%d pathogenic, %d benign)"
          % (l["n"], l["pathogenic"], l["benign"]))
    print("\n            pathogenic  benign")
    print("  firing    %10d  %6d" % (t["firing_pathogenic"], t["firing_benign"]))
    print("  silent    %10d  %6d" % (t["silent_pathogenic"], t["silent_benign"]))
    f = result["fisher"]
    print("\nFisher exact: OR %s  p = %.3g"
          % ("infinite" if f["odds_ratio_is_infinite"] else f["odds_ratio"], f["p_value"]))
    s = result["pathogenic_but_silent"]
    print("pathogenic but silent: %d variants across %d of %d systems"
          % (s["n"], s["systems_spanned"], s["systems_total"]))
    print("  %s" % ", ".join(s["variants"]))

    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print("\nwrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
