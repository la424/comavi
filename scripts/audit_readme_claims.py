#!/usr/bin/env python
"""Gate the README's headline numbers against generated reference outputs.

WHY THIS EXISTS
===============
The README is the most-read file in the repository and, until this script, no
audit read it at all. Every headline the manuscript reports is restated there --
both agreement denominators, mechanism-consistency, the four-way per-axis
decomposition, the tier gradient and its rank correlation -- so a value could be
corrected in the paper and left stale in the README with nothing to catch it.
That is the same coverage hole that had left the Abstract ungated, on a file
that more people will read than the paper.

Every expectation below is READ FROM a committed reference output, never typed
here. The literals are matched with ``audit_match.check_literal``, so a value
stated to greater precision than the data supports fails, and a value that is
correct in one mention but stale in another is reported rather than masked by
the good occurrence.

Exit 0 = every gated README claim matches generated data.
"""
from __future__ import annotations

import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from audit_match import check_literal, weak_needle  # noqa: E402
import apply_concordance_v5 as ac  # noqa: E402

README = REPO / "README.md"
LEDGER = REPO / "reference_outputs" / "COMAVI_numbers_ledger.json"
CANONICAL = REPO / "reference_outputs" / "scored_61var_canonical.csv"
ISDS_SUMMARY = REPO / "reference_outputs" / "isds_v1" / "ISDS_v1_summary.json"


def results_block(text: str) -> str:
    """The '## Benchmark results' section only.

    Scoping matters: '61' or '0.72' appear all over a README for unrelated
    reasons, and a gate that matched them anywhere would be satisfied by
    coincidence rather than by the claim it is meant to police.
    """
    i = text.index("## Benchmark results")
    j = text.find("\n## ", i + 5)
    seg = text[i:j if j > 0 else len(text)]
    # Flatten whitespace: markdown hard-wraps prose, so a phrase needle can
    # straddle a line break and a raw containment test would miss it.
    return " ".join(seg.split())


def main() -> int:
    if not README.exists():
        print("README.md absent -- nothing to audit", file=sys.stderr)
        return 0

    led = json.loads(LEDGER.read_text())
    isds = json.loads(ISDS_SUMMARY.read_text())
    blk = results_block(README.read_text())

    canonical = pd.read_csv(
        CANONICAL,
        low_memory=False,
    )

    partners = [
        partner
        for partner in ac.discover_partners(canonical)
        if partner
        and "_ci95_" not in partner
        and "_distinguishable_" not in partner
    ]

    axis_names = (
        "tier",
        "monomer",
        "fold",
        "binding",
    )

    def tally(mask):
        totals = {
            axis: [0, 0]
            for axis in axis_names
        }

        for _, row in canonical.loc[mask].iterrows():
            result = ac.structural_agreement_by_axis(
                row,
                partners,
                2.5,
                2.5,
                2.5,
            )

            for axis, (numerator, denominator) in result.items():
                totals[axis][0] += int(numerator)
                totals[axis][1] += int(denominator)

        return totals

    all_mask = pd.Series(
        True,
        index=canonical.index,
    )

    primary_mask = canonical[
        "mech_consistency_t25"
    ].isin(
        {
            "consistent",
            "partial",
            "inconsistent",
        }
    )

    axes = tally(all_mask)
    primary_axes = tally(primary_mask)

    derived_total = [
        sum(pair[0] for pair in axes.values()),
        sum(pair[1] for pair in axes.values()),
    ]

    primary_total = [
        sum(pair[0] for pair in primary_axes.values()),
        sum(pair[1] for pair in primary_axes.values()),
    ]

    ledger_total = [
        int(value)
        for value in led["SA_total"]
    ]

    ledger_axes = {
        axis: [
            int(led["SA_by_axis"][axis][0]),
            int(led["SA_by_axis"][axis][1]),
        ]
        for axis in axis_names
    }

    for axis in axis_names:
        entry = led["SA_by_axis"][axis]

        if len(entry) < 3:
            raise SystemExit(
                f"COMAVI_numbers_ledger.json SA_by_axis[{axis!r}] "
                "lacks its stored rate."
            )

        expected_rate = round(
            axes[axis][0] / axes[axis][1],
            4,
        )

        observed_rate = round(
            float(entry[2]),
            4,
        )

        if observed_rate != expected_rate:
            raise SystemExit(
                "COMAVI_numbers_ledger.json axis rate is stale: "
                f"axis={axis}, observed={observed_rate}, "
                f"derived={expected_rate}"
            )

    if ledger_total != derived_total:
        raise SystemExit(
            "COMAVI_numbers_ledger.json SA_total is stale: "
            f"ledger={ledger_total}, derived={derived_total}"
        )

    if ledger_axes != axes:
        raise SystemExit(
            "COMAVI_numbers_ledger.json SA_by_axis is stale: "
            f"ledger={ledger_axes}, derived={axes}"
        )

    if primary_total != [99, 132]:
        raise SystemExit(
            "Unexpected primary structural-agreement total: "
            f"{primary_total}"
        )

    sa_ok, sa_all = derived_total

    derived_sa = round(
        sa_ok / sa_all,
        4,
    )

    observed_sa = round(
        float(led["SA"]),
        4,
    )

    if observed_sa != derived_sa:
        raise SystemExit(
            "COMAVI_numbers_ledger.json scalar SA is stale: "
            f"observed={observed_sa}, derived={derived_sa}"
        )

    grad = led["tier_gradient"]
    rho, p = led["tier_spearman"]
    mc_n = led["MC_n"]

    observed = isds["observed"]
    positive_count = int(
        isds["population"]["positive"]
    )

    top_by_k = {
        int(row["k"]): row
        for row in isds["top_k"]
    }

    # The graded-population denominator is the all-rows denominator minus the
    # axes contributed by variants excluded as ungraded-by-construction. The
    # README states both conventions; derive the graded one rather than typing
    # it, so the pair cannot silently disagree with the ledger.
    sa_graded_den = primary_total[1]

    required = [
        # --- headline agreement, both denominator conventions ---
        (f"{sa_ok}/{sa_graded_den}", "structural agreement, graded population"),
        (f"{sa_ok}/{sa_all} = {sa_ok / sa_all:.3f}", "structural agreement, all rows"),
        (f"graded n={mc_n}", "mechanism-consistency population"),
        (f"({led['MC_tier_ablated']:.4f}, graded n={mc_n})", "mechanism-consistency value"),

        # --- per-axis decomposition: all four must reconcile to the headline ---
        (f"tier {axes['tier'][0]}/{axes['tier'][1]}", "per-axis tier"),
        (f"monomer-fold \u0394\u0394G {axes['monomer'][0]}/{axes['monomer'][1]}",
         "per-axis monomer-fold"),
        (f"complex-fold \u0394\u0394G {axes['fold'][0]}/{axes['fold'][1]}",
         "per-axis complex-fold"),
        (f"binding \u0394\u0394G {axes['binding'][0]}/{axes['binding'][1]}",
         "per-axis binding"),
        (f"denominators sum to {sa_all}", "per-axis denominators reconcile"),

        # --- fixed cohort-independent ISDS-v1 outputs ---
        (
            f"ISDS-v1 ROC AUC "
            f"{observed['isds_roc_auc']:.3f}",
            "ISDS-v1 ROC AUC",
        ),
        (
            f"average precision "
            f"{observed['isds_average_precision']:.3f}",
            "ISDS-v1 average precision",
        ),
        (
            f"Top 10: "
            f"{int(top_by_k[10]['structural_mechanisms_in_top_k'])} "
            f"structural-mechanism variants",
            "ISDS-v1 top-10 utility",
        ),
        (
            f"Top 20: "
            f"{int(top_by_k[20]['structural_mechanisms_in_top_k'])}"
            f"/{positive_count}",
            "ISDS-v1 top-20 recovery",
        ),
        (
            f"Top 25: "
            f"{int(top_by_k[25]['structural_mechanisms_in_top_k'])}"
            f"/{positive_count}",
            "ISDS-v1 top-25 recovery",
        ),

        # --- tier gradient: the monotonic pathogenicity claim ---
        (f"Tier 1 {grad['Tier 1'][2] * 100:.0f}%", "tier 1 pathogenic rate"),
        (f"Tier 2 {grad['Tier 2'][2] * 100:.0f}%", "tier 2 pathogenic rate"),
        (f"Tier 3 {grad['Tier 3'][2] * 100:.0f}%", "tier 3 pathogenic rate"),
        (f"Tier 4 {grad['Tier 4'][2] * 100:.0f}%", "tier 4 pathogenic rate"),
        (f"\u03c1 = \u2212{abs(rho):.2f}", "tier gradient rank correlation"),
        (f"p = {p}", "tier gradient p-value"),
    ]

    # The superseded 44-variant gradient must never reappear. The ledger itself
    # records these as SUPERSEDED; this turns that note into a check.
    forbidden = [
        ("12/12", "superseded 44-set Tier 1 count"),
        ("13/16", "superseded 44-set Tier 2 count"),
        ("OR = 6.48", "superseded underpowered odds ratio"),
        ("OR 6.48", "superseded underpowered odds ratio"),
    ]

    fails = []
    for needle, label in required:
        weak = weak_needle(needle)
        if weak:
            fails.append(f"UNGATEABLE [{label}]: {weak}")
            continue
        err = check_literal(blk, needle)
        if err:
            fails.append(f"MISSING in README results block [{label}]: {err}")

    for needle, label in forbidden:
        if needle in blk:
            fails.append(f"FORBIDDEN in README [{label}]: {needle!r} is superseded")

    print(f"\n{len(required)} required literals, {len(forbidden)} forbidden, "
          f"{len(fails)} failure(s)")
    for f in fails:
        print(f"  {f}", file=sys.stderr)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
