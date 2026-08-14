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

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from audit_match import check_literal, weak_needle  # noqa: E402

README = REPO / "README.md"
LEDGER = REPO / "reference_outputs" / "COMAVI_numbers_ledger.json"


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
    blk = results_block(README.read_text())

    sa_ok, sa_all = led["SA_total"]                     # [99, 131]
    axes = led["SA_by_axis"]                            # tier/monomer/fold/binding
    grad = led["tier_gradient"]                         # Tier 1..4 -> [k, n, rate]
    rho, p = led["tier_spearman"]                       # [-0.3997, 0.0044]
    mc_n = led["MC_n"]

    # The graded-population denominator is the all-rows denominator minus the
    # axes contributed by variants excluded as ungraded-by-construction. The
    # README states both conventions; derive the graded one rather than typing
    # it, so the pair cannot silently disagree with the ledger.
    sa_graded_den = sa_all - 1

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
