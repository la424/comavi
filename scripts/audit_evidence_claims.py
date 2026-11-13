#!/usr/bin/env python
"""Assert every evidence-ledger number in the manuscript matches the generated data.

The §2.2 and §3.9 prose asserts ~30 literals about the evidence ledger. Each is
checked here against the two JSON summaries so a re-curation cannot silently
desynchronise the paper from its data. Prose is wrong until the data says
otherwise.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from audit_match import check_literal, weak_needle  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parent.parent
MS = REPO / "docs" / "COMAVI_manuscript_current.txt"
LED = REPO / "reference_outputs" / "COMAVI_evidence_ledger_summary.json"
STRAT = REPO / "reference_outputs" / "COMAVI_evidence_stratified_agreement.json"
LEDGER = REPO / "reference_outputs" / "COMAVI_numbers_ledger.json"
TIER = REPO / "reference_outputs" / "COMAVI_tier_construction.json"
LEDGER_CSV = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"

E1 = "E1_quantitative_energetic"
E2 = "E2_quantitative_functional"
E3 = "E3_qualitative_experimental"
E4 = "E4_population_frequency"
E5 = "E5_inferred_no_axis_assay"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--allow-missing-manuscript",
        action="store_true",
        help="Skip instead of failing when the manuscript text is absent. For "
             "public clones that withhold the manuscript ONLY. Never pass this "
             "in the working repository: a silent skip here is exactly how the "
             "prose drifted from the data undetected.",
    )
    args = ap.parse_args()

    if not MS.exists():
        if args.allow_missing_manuscript:
            print(f"[SKIP] {MS.relative_to(REPO)} absent; skip explicitly requested.")
            return 0
        print(f"[FAIL] {MS.relative_to(REPO)} is not present.", file=sys.stderr)
        print("       Regenerate it with scripts/extract_manuscript_text.py, or pass",
              file=sys.stderr)
        print("       --allow-missing-manuscript if this is a manuscript-withheld clone.",
              file=sys.stderr)
        return 1
    text = MS.read_text()
    led = json.loads(LED.read_text())
    st = json.loads(STRAT.read_text())
    nl = json.loads(LEDGER.read_text())
    tier = json.loads(TIER.read_text())["screen_full_tier"]
    # v11: the SI states the benchmark's 97 expectations, not the whole 99-row
    # ledger (two rows belong to excluded BRCT candidates), so the evidence
    # composition is counted on the benchmark subset of the committed ledger.
    sys.path.insert(0, str(REPO / "scripts"))
    import benchmark_population as bp
    _can = pd.read_csv(REPO / "reference_outputs" / "scored_61var_canonical.csv", low_memory=False)
    _bk = set(zip(bp.benchmark(_can).system, bp.benchmark(_can).variant))
    _led = pd.read_csv(LEDGER_CSV)
    _lb = _led[[k in _bk for k in zip(_led.system, _led.variant)]]
    assert len(_lb) == json.loads(STRAT.read_text())["accounting"]["n_committed_axes"], len(_lb)
    typ = _lb["evidence_type"].value_counts().to_dict()
    dirn = _lb["evidence_directness"].value_counts().to_dict()
    acc, crude, strat = st["accounting"], st["crude"], st["stratified"]
    by_type = st["by_evidence_type"]
    comp = st["token_composition_by_grade"]
    _fcb = _lb[_lb["axis"] == "fold_complex"]
    fold_inferred = int((_fcb["evidence_type"] == E5).sum())
    fold_committed = int(len(_fcb))
    # Directness-inferred share of the complex-fold axis. Not in the summary
    # JSON, so it is read from the committed ledger. Note S2 states BOTH this
    # (14/25, directness) and fold_inferred (11, evidence type); they are
    # different quantities and conflating them is how that sentence read
    # ambiguously through v7.7.
    fold_inferred_directness = int((_fcb["evidence_directness"] == "inferred").sum())

    _amc = pd.read_csv(REPO / "reference_outputs" / "isds_v1" /
                       "ISDS_v1_alphamissense_pathogenicity_comparison.csv").set_index("descriptor")
    am_auc = {k: "%.3f" % _amc.loc[k, "pathogenicity_auc"] for k in ("AM pathogenicity", "isds_v1")}
    _p4 = "%.4f" % tier["within_system_permutation_p"]
    assert _p4[-1] != "5", "stored permutation p ends in 5; three-place rounding would be ambiguous"
    perm3 = "%.3f" % float(_p4)

    def _mh(per_token):
        num = den = 0.0
        for k in ("destab", "neutral"):
            s_ = per_token[k]
            a, b = s_["quant_ok"], s_["quant_n"] - s_["quant_ok"]
            c, d = s_["soft_ok"], s_["soft_n"] - s_["soft_ok"]
            n_ = a + b + c + d
            num += a * d / n_
            den += b * c / n_
        return num / den
    mh_or = _mh(strat["per_token"])
    assert abs(mh_or - strat["OR_mantel_haenszel"]) < 5e-4, (mh_or, strat["OR_mantel_haenszel"])

    checks = [
        # Needles are built from GENERATED values and asserted to appear in the
        # manuscript text. They were keyed to a §2.2 prose enumeration that no
        # longer exists (the evidence composition moved into Table S2 and Note
        # S2), which is why this gate reported nothing for several revisions.
        # Every needle below is verified present against the adjudicated draft.
        # ---- Note S2 ledger scope and composition
        # v11: S1 Text states the benchmark's expectations (97 on the 56 graded
        # variants), not the whole ledger's, so the count is read from the
        # benchmark accounting of the stratified record.
        (f"variants carry {acc['n_committed_axes']} committed energetic-axis expectations",
         "ledger scope"),
        (f"all {typ[E4]} E4 axes and all {typ[E5]} E5 axes are neutral",
         "E4/E5 neutrality"),
        (f"{fold_inferred_directness} of its {fold_committed} committed axes "
         f"rest on inferred", "fold-complex inferred directness"),
        (f"and {fold_inferred} are E5 inferred axis states",
         "fold-complex E5 count"),
        # ---- Table S2 evidence composition (tab-separated; the extractor
        # renders Word tables with tabs, so pipe-delimited needles never match)
        (f"E1 quantitative energetic\t{typ[E1]}\t", "Table S2 E1 axes"),
        (f"E2 quantitative functional\t{typ[E2]}\t", "Table S2 E2 axes"),
        (f"E3 qualitative experimental\t{typ[E3]}\t", "Table S2 E3 axes"),
        (f"E4 population frequency\t{typ[E4]}\t", "Table S2 E4 axes"),
        (f"E5 inferred (no axis assay)\t{typ[E5]}\t", "Table S2 E5 axes"),
        # ---- Note S2 crude contrast
        (f"Of the {acc['n_committed_axes']} expectations, {acc['n_ddg_axes_gradeable']} can be "
         f"compared with a prediction", "gradeable energetic axes"),
        (f"{crude['quant_ok']}/{crude['quant_n']} = {crude['quant_rate']:.3f} "
         f"for E1\u2013E2 evidence", "crude quant"),
        (f"{crude['soft_ok']}/{crude['soft_n']} = {crude['soft_rate']:.3f} "
         f"for E3\u2013E5 evidence", "crude soft"),
        (f"(Fisher p = {crude['p']:.3f})", "crude Fisher p"),
        # ---- Note S2 stratified contrast. The null result is the CLAIM here,
        # so both strata and the Mantel-Haenszel estimate must be gated: a
        # stale stratum that still reads null would pass an eyeball check.
        (f"agreement was {strat['per_token']['destab']['quant_ok']}/"
         f"{strat['per_token']['destab']['quant_n']} for E1\u2013E2 and "
         f"{strat['per_token']['destab']['soft_ok']}/"
         f"{strat['per_token']['destab']['soft_n']} for E3\u2013E5", "destab stratum"),
        (f"the corresponding rates were "
         f"{strat['per_token']['neutral']['quant_ok']}/"
         f"{strat['per_token']['neutral']['quant_n']} and "
         f"{strat['per_token']['neutral']['soft_ok']}/"
         f"{strat['per_token']['neutral']['soft_n']}", "neutral stratum"),
        # The stored odds ratio is rounded to 3 places (0.795), and formatting
        # that to 2 rounds twice (0.80); the exact estimate from the stratum
        # counts is 0.7949, so the needle is built from the counts.
        (f"odds ratio across expected-state strata was {mh_or:.2f} "
         f"(continuity-corrected test p = {strat['p']:.2f})", "Mantel-Haenszel"),
        (f"Among {by_type[E1]['n']} gradeable E1 energetic axes, "
         f"{by_type[E1]['ok']} agreed ({by_type[E1]['rate']:.3f})",
         "E1 agreement"),
        # ---- §3.14 the AlphaMissense comparison must state BOTH AUCs.
        # The accuracy gap is the premise of the orthogonality argument, not a
        # concession to it: quoting AM's 0.903 while dropping the tier's 0.769
        # turns an honest comparison into an unfalsifiable claim of
        # complementarity. A v30-lineage draft did exactly that and passed every
        # check, because these two literals were gated only in the deck audit.
        # v11: the comparison is on the 39 labelled prioritization variants
        # (ISDS_v1_alphamissense_pathogenicity_comparison.csv), and the score
        # set against AlphaMissense is the priority score. The legacy
        # AM_auc / tier_auc_on_AM_set keys of the numbers ledger describe the
        # retired 47-variant set and are no longer read.
        (am_auc["AM pathogenicity"], "AlphaMissense AUC"),
        (am_auc["isds_v1"], "priority-score AUC on the AlphaMissense set (the gap)"),
        # ---- §3.12 tier screen: the system-level conditioning is what makes the
        # perfect-recall screen more than a small-sample accident. Report the
        # EXACT enumerated permutation p from verify_tier_construction.py, never
        # the 20k Monte-Carlo estimate in analyze_tier_energy_gating.py.
        # reported to three places in v11; the record stores four, so refuse a tie
        (perm3, "exact within-system permutation p"),
        (f"{tier['sensitivity_ci'][0]}", "screen sensitivity CI lower bound"),
        (f"{tier['specificity_ci'][0]}", "screen specificity CI lower bound"),
    ]

    # Boundary-aware matching (see scripts/audit_match.py): a plain `needle in
    # text` test passes when prose states a MORE precise value than the data
    # supports, because "0.903" is contained in "0.9039". Mutation-testing found
    # every digit-ending gate here non-binding under plain containment.
    fails = []
    for needle, label in checks:
        why = weak_needle(needle)
        if why is not None:
            print(f"UNGATEABLE [{label}]: {why}", file=sys.stderr)
            fails.append(label)
            continue
        problem = check_literal(text, needle)
        if problem is not None:
            print(f"MISSING [{label}]: {problem}", file=sys.stderr)
            fails.append(label)

    # Forbidden: the retired prose that claimed evidence grading without data.
    forbidden = [
        ("graded by type: **quantitative**", "pre-ledger evidence-grading claim"),
        ("structural-inference** (contact-level", "pre-ledger evidence tiers"),
    ]
    for needle, label in forbidden:
        if needle in text:
            print(f"FORBIDDEN [{label}] still present", file=sys.stderr)
            fails.append(label)

    print(f"evidence-claim audit: {len(checks)-len([f for f in fails if f not in dict(forbidden).values()])}"
          f"/{len(checks)} literals verified")
    if fails:
        print(f"FAIL ({len(fails)})", file=sys.stderr)
        return 1
    print("PASS — every evidence-ledger literal in the manuscript matches generated data")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
