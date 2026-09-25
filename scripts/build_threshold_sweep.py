#!/usr/bin/env python3
"""Every reported metric as a function of the decision threshold.

WHY THIS EXISTS
---------------
COMAVI's outputs are continuous: three per-axis energies, a structural-context
tier, and a priority score. None of them needs a threshold. The threshold
exists only to convert continuous energies into the discrete mechanism CALLS
that can be graded against categorical literature labels, so it is an
instrument of this evaluation rather than a parameter of the pipeline.

That distinction has a consequence the paper has to face: every graded number
is a function of one evaluation choice. Reporting a single operating point
invites the reader to treat it as the method's accuracy. This generator emits
the whole curve so the reference point can be read as one column of a table
rather than as the answer.

WHAT IT EMITS, AT EACH OF THE FIVE THRESHOLDS
---------------------------------------------
  detection          of curated structural mechanisms, how many have ANY axis
                     firing
  attribution        of those, how many have the EXPECTED axis firing
  correct_rejection  of curated no-lesion variants, how many have no axis
                     firing
  whole_variant      the partial-credit mechanism-pattern score, from the
                     canonical's stored per-threshold grade columns
  per_axis           directional agreement for monomer fold, complex-context
                     fold, binding and tier SEPARATELY, delegated to
                     apply_concordance_v5.structural_agreement_by_axis so the
                     decomposition cannot drift from the headline it sums to
  structural_agreement  the pooled four-way total

WHY NOT A FINER OR WIDER GRID
-----------------------------
The five thresholds are the committed set in apply_concordance_v5.
THRESHOLD_SPECS, unchanged since the repository's first commit. Four are scalar
(1.0, 1.5, 2.0, 2.5 kcal/mol); the fifth is per-axis (2.9 monomer, 2.9
complex-context fold, 3.5 binding). Intermediate or higher points are
computable -- every quantity here is a reduction over stored energies -- but
they are NOT in the committed set, and adding them would mean grading the
benchmark at operating points chosen after the ground truth was fixed. The
sweep therefore stops where the prespecified set stops, and `--extended`
reports the additional points separately, clearly labelled post hoc, so the
top-of-range behaviour can be inspected without being promoted into the
primary table.

Usage:
    PYTHONPATH=. python scripts/build_threshold_sweep.py
    PYTHONPATH=. python scripts/build_threshold_sweep.py --check
    PYTHONPATH=. python scripts/build_threshold_sweep.py --extended
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
# Sibling generators are imported by module name, so scripts/ must be on the
# path regardless of how PYTHONPATH is set. CI sets it to the repo root only,
# under which the bare imports below raise ModuleNotFoundError -- this gate
# passed locally and would have failed in CI without this line.
sys.path.insert(0, str(REPO / "scripts"))
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT_JSON = REPO / "reference_outputs" / "COMAVI_threshold_sweep.json"
OUT_CSV = REPO / "reference_outputs" / "COMAVI_threshold_sweep.csv"

# The post hoc points reported under --extended get their OWN committed record
# rather than extra rows in the primary one. Two reasons. The primary record's
# --check compares stored rows against a build with the prespecified set only,
# so extended rows there would make the gate fail against itself. And the
# supplement publishes the post hoc columns: before this file existed, those
# values appeared in S2 Table with nothing in the repository to check them
# against, which is the same defect class as a hardcoded literal. Anything
# published from --extended must be reproducible from a committed file.
OUT_EXT_JSON = REPO / "reference_outputs" / "COMAVI_threshold_sweep_extended.json"
OUT_EXT_CSV = REPO / "reference_outputs" / "COMAVI_threshold_sweep_extended.csv"

GRADE = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
N_BOOT = 2000
SEED = 20260923

# Post hoc points, reported only under --extended and never in the primary
# table. Named here so the set is auditable rather than ad hoc.
EXTENDED = [("t30", 3.0), ("t35", 3.5), ("t40", 4.0)]


def _specs():
    """The committed threshold set, as (tag, spec) with spec scalar or per-axis."""
    import apply_concordance_v5 as ac

    return list(ac.THRESHOLD_SPECS)


def _per_axis(spec):
    """Expand a spec into (monomer, fold, binding) thresholds."""
    if isinstance(spec, dict):
        return spec["monomer"], spec["fold"], spec["binding"]
    return float(spec), float(spec), float(spec)


def sweep_row(canon, partners, tag, spec, ac, dva):
    """One threshold's worth of every reported quantity."""
    mono, fold, bind = _per_axis(spec)
    graded = canon[canon.mech_consistency_t25.isin(GRADE)].copy()

    # --- detection / attribution / correct rejection -----------------------
    # Delegated to the shipped generator's own build(), parameterised by
    # threshold, so the sweep is the same computation at five points rather
    # than a second implementation that happens to agree at one.
    built = dva.build(canon, threshold=spec)

    # --- whole-variant score ----------------------------------------------
    col = "mech_consistency_%s" % tag
    if col in canon.columns:
        g = canon[canon[col].isin(GRADE)]
        wv_k = float(g[col].map(GRADE).sum())
        wv_n = int(len(g))
    else:
        wv_k, wv_n = float("nan"), 0

    # --- per-axis directional agreement -----------------------------------
    per = {k: [0, 0] for k in ("tier", "monomer", "fold", "binding")}
    for _, row in graded.iterrows():
        d = ac.structural_agreement_by_axis(row, partners, mono, fold, bind)
        for k, (n, dd) in d.items():
            per[k][0] += n
            per[k][1] += dd

    return {
        "threshold_tag": tag,
        "threshold_kcal_mol": (
            "%.1f" % mono if mono == fold == bind
            else "%.1f/%.1f/%.1f" % (mono, fold, bind)),
        "per_axis_spec": mono == fold == bind,
        "detection_k": built["detection"]["k"],
        "detection_n": built["detection"]["n"],
        "detection": round(built["detection"]["k"] / built["detection"]["n"], 4),
        "attribution_k": built["attribution"]["k"],
        "attribution_n": built["attribution"]["n"],
        "attribution": round(built["attribution"]["k"] / built["attribution"]["n"], 4),
        "correct_rejection_k": built["correct_rejection"]["k"],
        "correct_rejection_n": built["correct_rejection"]["n"],
        "correct_rejection": round(
            built["correct_rejection"]["k"] / built["correct_rejection"]["n"], 4),
        "whole_variant_k": round(wv_k, 1),
        "whole_variant_n": wv_n,
        "whole_variant": round(wv_k / wv_n, 4) if wv_n else float("nan"),
        "axis_monomer_k": per["monomer"][0], "axis_monomer_n": per["monomer"][1],
        "axis_monomer": round(per["monomer"][0] / per["monomer"][1], 4) if per["monomer"][1] else float("nan"),
        "axis_fold_k": per["fold"][0], "axis_fold_n": per["fold"][1],
        "axis_fold": round(per["fold"][0] / per["fold"][1], 4) if per["fold"][1] else float("nan"),
        "axis_binding_k": per["binding"][0], "axis_binding_n": per["binding"][1],
        "axis_binding": round(per["binding"][0] / per["binding"][1], 4) if per["binding"][1] else float("nan"),
        "axis_tier_k": per["tier"][0], "axis_tier_n": per["tier"][1],
        "axis_tier": round(per["tier"][0] / per["tier"][1], 4) if per["tier"][1] else float("nan"),
        "structural_agreement_k": sum(v[0] for v in per.values()),
        "structural_agreement_n": sum(v[1] for v in per.values()),
        "structural_agreement": round(
            sum(v[0] for v in per.values()) / sum(v[1] for v in per.values()), 4),
    }


def cluster_band(canon, tag, n_boot=N_BOOT, seed=SEED):
    """System-level bootstrap on the whole-variant score at one threshold.

    Variants inside a system share a structure and a curation pass, so the
    resampling unit is the system, not the variant.
    """
    col = "mech_consistency_%s" % tag
    if col not in canon.columns:
        return None
    g = canon[canon[col].isin(GRADE)][["system", col]].copy()
    g["s"] = g[col].map(GRADE)
    systems = g.system.unique()
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        pick = rng.choice(systems, size=len(systems), replace=True)
        vals = np.concatenate([g.loc[g.system == s, "s"].to_numpy() for s in pick])
        draws.append(vals.mean())
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return {"lo": round(float(lo), 4), "hi": round(float(hi), 4),
            "sd": round(float(np.std(draws)), 4), "n_boot": n_boot,
            "n_systems": int(len(systems))}


def build(extended=False):
    import apply_concordance_v5 as ac
    import build_detection_vs_attribution as dva

    canon = pd.read_csv(CANON, low_memory=False)
    partners = ac.discover_partners(canon)

    specs = _specs()
    if extended:
        specs = specs + EXTENDED

    rows = [sweep_row(canon, partners, tag, spec, ac, dva) for tag, spec in specs]
    table = pd.DataFrame(rows)
    bands = {tag: cluster_band(canon, tag) for tag, _ in _specs()}

    committed = [t for t, _ in _specs()]
    result = {
        "committed_thresholds": committed,
        "extended_thresholds_post_hoc": [t for t, _ in EXTENDED] if extended else [],
        "provenance": (
            "THRESHOLD_SPECS has been unchanged since the repository's first "
            "commit; the evidence ledger that defines the graded ground truth "
            "was first committed two months later. The threshold set was "
            "therefore fixed before the ground truth existed in committed form."),
        "rows": rows,
        "cluster_bands": bands,
    }
    return result, table


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="verify the committed record matches a fresh build")
    ap.add_argument("--extended", action="store_true",
                    help="also report post hoc thresholds above the committed set")
    args = ap.parse_args()

    result, table = build(extended=args.extended)
    ref = next(r for r in result["rows"] if r["threshold_tag"] == "t25")

    # The gate that matters: the reference column must reproduce every value
    # the manuscript currently publishes. If it does not, this sweep is a
    # second implementation rather than the same computation at five points.
    import build_detection_vs_attribution as dva
    pub = json.loads(
        (REPO / "reference_outputs" / "COMAVI_detection_vs_attribution.json").read_text())
    for key in ("detection", "attribution", "correct_rejection"):
        got = (ref["%s_k" % key], ref["%s_n" % key])
        want = (pub[key]["k"], pub[key]["n"])
        assert got == want, (
            "reference column disagrees with the published %s: sweep %s, "
            "published %s" % (key, got, want))
    assert ref["whole_variant_k"] == 41.0 and ref["whole_variant_n"] == 57, (
        "reference whole-variant score moved: %s/%s" % (
            ref["whole_variant_k"], ref["whole_variant_n"]))

    if args.check:
        target = OUT_EXT_JSON if args.extended else OUT_JSON
        if not target.exists():
            print("FAIL: %s missing" % target.name)
            return 1
        stored = json.loads(target.read_text())
        # The post hoc rows carry no whole-variant score, so that field is NaN.
        # A plain != comparison is always True in the presence of NaN, which
        # would make this gate fail against a record it had just written.
        def same(a, b):
            if isinstance(a, float) and isinstance(b, float):
                if a != a and b != b:
                    return True
                return abs(a - b) < 1e-12
            if isinstance(a, dict) and isinstance(b, dict):
                return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
            if isinstance(a, list) and isinstance(b, list):
                return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
            return a == b

        if not same(stored.get("rows"), result["rows"]):
            print("FAIL: stored sweep differs from a fresh build (%s)" % target.name)
            return 1
        if not args.extended and not OUT_EXT_JSON.exists():
            print("FAIL: %s missing; the supplement publishes the post hoc "
                  "columns and they must be reproducible" % OUT_EXT_JSON.name)
            return 1
        print("PASS: threshold sweep reproduces, and the reference column "
              "matches every published value")
        return 0

    cols = ["threshold_kcal_mol", "detection", "attribution", "correct_rejection",
            "whole_variant", "axis_monomer", "axis_fold", "axis_binding",
            "structural_agreement"]
    print("Every reported metric as a function of decision threshold\n")
    print(table[cols].to_string(index=False))
    print("\nwhole-variant score with system-level bootstrap bands:")
    for tag, b in result["cluster_bands"].items():
        r = next(x for x in result["rows"] if x["threshold_tag"] == tag)
        print("  %-6s %-12s %.4f  (95%% CI %.4f-%.4f, sd %.4f)"
              % (tag, r["threshold_kcal_mol"], r["whole_variant"],
                 b["lo"], b["hi"], b["sd"]))

    peak = max(result["rows"], key=lambda r: r["whole_variant"])
    print("\n  argmax of the whole-variant score: %s (%s kcal/mol) at %.4f"
          % (peak["threshold_tag"], peak["threshold_kcal_mol"], peak["whole_variant"]))

    out_json, out_csv = (OUT_EXT_JSON, OUT_EXT_CSV) if args.extended else (OUT_JSON, OUT_CSV)
    out_json.write_text(json.dumps(result, indent=1) + "\n")
    table.to_csv(out_csv, index=False)
    print("\nwrote %s" % out_json.relative_to(REPO))
    print("wrote %s" % out_csv.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
