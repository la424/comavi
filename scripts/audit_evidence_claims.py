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
    typ, dirn = led["by_evidence_type"], led["by_directness"]
    acc, crude, strat = st["accounting"], st["crude"], st["stratified"]
    by_type = st["by_evidence_type"]
    comp = st["token_composition_by_grade"]
    fold_inferred = led["type_by_axis"]["fold_complex"][E5]
    fold_committed = sum(led["type_by_axis"]["fold_complex"].values())
    # Directness-inferred share of the complex-fold axis. Not in the summary
    # JSON, so it is read from the committed ledger. Note S2 states BOTH this
    # (14/25, directness) and fold_inferred (11, evidence type); they are
    # different quantities and conflating them is how that sentence read
    # ambiguously through v7.7.
    _fc = pd.read_csv(LEDGER_CSV)
    _fc = _fc[_fc["axis"] == "fold_complex"]
    fold_inferred_directness = int((_fc["evidence_directness"] == "inferred").sum())

    checks = [
        # Needles are built from GENERATED values and asserted to appear in the
        # manuscript text. They were keyed to a §2.2 prose enumeration that no
        # longer exists (the evidence composition moved into Table S2 and Note
        # S2), which is why this gate reported nothing for several revisions.
        # Every needle below is verified present against the adjudicated draft.
        # ---- Note S2 ledger scope and composition
        (f"contains {led['n_committed_axes']} committed energetic-axis expectations",
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
        (f"E5 inferred\t{typ[E5]}\t", "Table S2 E5 axes"),
        # ---- Note S2 crude contrast
        (f"all-row {acc['n_ddg_axes_gradeable']}-axis convention",
         "gradeable energetic axes"),
        (f"{crude['quant_ok']}/{crude['quant_n']} = {crude['quant_rate']:.3f} "
         f"for E1-E2 evidence", "crude quant"),
        (f"{crude['soft_ok']}/{crude['soft_n']} = {crude['soft_rate']:.3f} "
         f"for E3-E5 evidence", "crude soft"),
        (f"(Fisher p = {crude['p']:.3f})", "crude Fisher p"),
        # ---- Note S2 stratified contrast. The null result is the CLAIM here,
        # so both strata and the Mantel-Haenszel estimate must be gated: a
        # stale stratum that still reads null would pass an eyeball check.
        (f"agreement was {strat['per_token']['destab']['quant_ok']}/"
         f"{strat['per_token']['destab']['quant_n']} for E1-E2 and "
         f"{strat['per_token']['destab']['soft_ok']}/"
         f"{strat['per_token']['destab']['soft_n']} for E3-E5", "destab stratum"),
        (f"the corresponding rates were "
         f"{strat['per_token']['neutral']['quant_ok']}/"
         f"{strat['per_token']['neutral']['quant_n']} and "
         f"{strat['per_token']['neutral']['soft_ok']}/"
         f"{strat['per_token']['neutral']['soft_n']}", "neutral stratum"),
        (f"odds ratio {strat['OR_mantel_haenszel']:.2f} "
         f"(p = {strat['p']:.2f})", "Mantel-Haenszel"),
        (f"Among {by_type[E1]['n']} gradeable E1 energetic axes, "
         f"{by_type[E1]['ok']} agreed ({by_type[E1]['rate']:.3f})",
         "E1 agreement"),
        # ---- §3.14 the AlphaMissense comparison must state BOTH AUCs.
        # The accuracy gap is the premise of the orthogonality argument, not a
        # concession to it: quoting AM's 0.903 while dropping the tier's 0.769
        # turns an honest comparison into an unfalsifiable claim of
        # complementarity. A v30-lineage draft did exactly that and passed every
        # check, because these two literals were gated only in the deck audit.
        (f"{nl['AM_auc']}", "AlphaMissense AUC"),
        (f"{nl['tier_auc_on_AM_set']}", "tier AUC on the AM set (the gap)"),
        # ---- §3.12 tier screen: the system-level conditioning is what makes the
        # perfect-recall screen more than a small-sample accident. Report the
        # EXACT enumerated permutation p from verify_tier_construction.py, never
        # the 20k Monte-Carlo estimate in analyze_tier_energy_gating.py.
        (f"{tier['within_system_permutation_p']}", "exact within-system permutation p"),
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
