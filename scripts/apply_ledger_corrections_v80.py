#!/usr/bin/env python3
"""v8.0 -- Five axis-expectation corrections verified against primary sources.

WHY THIS EXISTS
---------------
An external re-curation (R1) proposed six changes to committed axis
expectations. Each was checked here against the primary source itself rather
than against the re-curation's summary of it. Five are applied. One is not.

APPLIED -- and what the source actually says
--------------------------------------------
kras_craf G12D  binding  neutral -> destab
kras_craf G12V  binding  neutral -> destab
    Hunter 2015 Table 2 ("Relative affinities for RAF kinase", nmol/L):
    WT 56 +/- 6, G12D 270 +/- 46, G12V 411 +/- 40. The text reports a
    "4.8-, 7.3-, and 6.2-fold decrease in affinity" for D, V and R.

    These two cells were changed to destab in v7.6, reverted to neutral in
    v7.9, and are now destab again. That is not churn. v7.9 reverted them
    because the row's own basis read "WT-like RAF1 binding", and a token
    contradicting its basis is exactly the defect v7.9 existed to fix. But the
    basis itself was wrong: it cited Hunter 2015, and Hunter 2015 measures a
    decrease. v7.9 corrected the token to agree with a basis that misread its
    own source. This correction fixes the basis and the token together.

    Methodological consequence, worth stating in Limitations: a gate that
    checks token against basis cannot catch a basis that misreads its source.
    verify_ledger_token_basis.py passed on these rows throughout.

mlh1_pms2 H718Y  monomer  (new row) destab
    Hinrichsen 2013: "the 4 MMR-proficient variants (K618A, E578G, V716M, and
    H718Y) moderately destabilize the MLH1 protein but are not causative for
    Lynch syndrome". Thermal stability "decreased more strongly in V716M,
    H718Y, and E578G"; pulse-chase half-life reduced to 64% on average.
    Structurally destabilizing and clinically neutral in one variant, from the
    source's own words.

smad4_smad3 I500V  binding  (new row) stab
    Lindsay 2025 Fig 4D-E: GST pull-down shows increased recovery of SMAD3
    with SMAD4 I500V, and TR-FRET shows a statistically significant increase
    in SMAD4-SMAD3 interaction. R361H and R361C were run as negative controls
    and show the opposite direction, so the assay resolves sign.

msh2_msh6 G674R  binding  neutral -> UNCOMMITTED (row removed)
    No accessible source establishes a binding direction for this exact
    substitution. The v7.15 basis reads "MSH6 binding preserved by co-IP" and
    cites Jia 2021, which measures MMR loss of function by 6-thioguanine
    selection and performs no co-immunoprecipitation. Ollila 2008 did perform
    co-IP on an MSH2 ATPase-domain variant, but on G674A -- a different
    substitution at the same position; that paper never mentions G674R.
    Removing the commitment is the conservative reading.

NOT APPLIED
-----------
msh2_msh6 C697F  binding  neutral -> destab   REJECTED, evidence insufficient
    R1 cites Lutzen 2008 (closed access, unavailable here). The existing
    neutral commitment is affirmatively supported by a source that could be
    read: Ollila 2006 Gastroenterology tested C697F by co-immunoprecipitation
    and reports "None of the studied MSH2 mutations destroyed the protein or
    abolished MSH2/MSH6 interaction"; C697F's summary row reads Normal
    expression, Normal interaction, Deficient repair. Codon 697 also lies
    outside both MSH6-interaction regions that paper identifies (378-625,
    875-934). Overturning a source that was read, on the strength of one that
    was not, is not a correction.

    Note that this rejection costs us: C697F was one of only two proposed
    changes that RAISED the mechanism-pattern score and the margin over the
    trivial baseline. It is held anyway. Revisit if Lutzen 2008 becomes
    available.

Ledger arithmetic: 98 - 1 (G674R) + 2 (H718Y monomer, I500V binding) = 99.

Usage:
    python scripts/apply_ledger_corrections_v80.py [--dry-run] [--check]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

# Import the derivation path and cell renderer rather than copying them, so
# this correction cannot drift from the ones that preceded it.
from apply_ledger_corrections_v76 import _cell, rederive  # noqa: E402

CANONICAL = ROOT / "reference_outputs" / "scored_61var_canonical.csv"
LEDGER = ROOT / "reference_outputs" / "COMAVI_evidence_ledger.csv"
SNAPDIR = ROOT / "reference_outputs"

HUNTER = ("Hunter 2015 Mol Cancer Res 13:1325 Table 2 "
          "(doi:10.1158/1541-7786.MCR-15-0203)")
HINRICHSEN = ("Hinrichsen 2013 Clin Cancer Res 19:2432 "
              "(doi:10.1158/1078-0432.CCR-12-3299)")
LINDSAY = ("Lindsay 2025 J Allergy Clin Immunol Fig 4D-E "
           "(doi:10.1016/j.jaci.2024.08.024)")

# system, variant, axis, from_token, to_token, type, directness, basis, citation
UPDATES = [
    ("kras_craf", "G12D", "binding", "neutral", "destab",
     "E1_quantitative_energetic", "direct",
     "RAF-RBD affinity decreased 4.8-fold (270 +/- 46 nmol/L against 56 +/- 6 "
     "for wild type); the source reports a decrease in affinity. Measured "
     "magnitude 0.93 kcal/mol, below the 2.5 kcal/mol reference threshold.",
     HUNTER),
    ("kras_craf", "G12V", "binding", "neutral", "destab",
     "E1_quantitative_energetic", "direct",
     "RAF-RBD affinity decreased 7.3-fold (411 +/- 40 nmol/L against 56 +/- 6 "
     "for wild type); the source reports a decrease in affinity. Measured "
     "magnitude 1.18 kcal/mol, below the 2.5 kcal/mol reference threshold.",
     HUNTER),
]

ADDITIONS = [
    ("mlh1_pms2", "H718Y", "monomer", "destab",
     "E2_quantitative_functional", "direct",
     "Source states the MMR-proficient variants including H718Y moderately "
     "destabilize the MLH1 protein; thermal stability decreased more strongly "
     "in V716M, H718Y and E578G, and pulse-chase half-life fell to 64% of "
     "wild type on average. Clinically classified neutral.",
     HINRICHSEN),
    ("smad4_smad3", "I500V", "binding", "stab",
     "E2_quantitative_functional", "direct",
     "GST pull-down shows increased recovery of SMAD3 with SMAD4 I500V, and "
     "TR-FRET shows a statistically significant increase in SMAD4-SMAD3 "
     "interaction signal. R361H and R361C negative controls show the opposite "
     "direction, so the assay resolves sign.",
     LINDSAY),
]

REMOVALS = [
    ("msh2_msh6", "G674R", "binding", "neutral",
     "no accessible source establishes a binding direction for this exact "
     "substitution; the prior basis attributed a co-IP result to a paper that "
     "performed no co-IP, and the co-IP that exists was run on G674A"),
]

COL = {"monomer": "expected_ddg_monomer",
       "fold_complex": "expected_ddg_fold_complex",
       "binding": "expected_ddg_binding"}

REJECTED = [("msh2_msh6", "C697F", "binding", "neutral",
             "Lutzen 2008 unavailable; Ollila 2006 co-IP reports interaction "
             "intact for this exact variant")]


def _match(led, system, variant, axis):
    return ((led["system"] == system) & (led["variant"] == variant)
            & (led["axis"] == axis))


# The ledger is NOT patched here. It is regenerated from the canonical by
# scripts/curate_evidence_types.py, whose ASSIGN table carries the corrected
# evidence_type, directness, basis and citation for each of these axes. That
# keeps one source of truth: this script moves the expectation cells in the
# canonical, and the generator derives the ledger from them. Run order is
#   apply_ledger_corrections_v80.py  ->  curate_evidence_types.py
# Earlier correction scripts (v7.6, v7.9) patched the ledger CSV directly;
# doing that here as well would leave the generator and the patch able to
# disagree, and would make a clean checkout non-reproducible because the patch
# asserts a from-token the generator has already moved.


def patch_canonical(dry_run: bool) -> pd.DataFrame:
    base = pd.read_csv(CANONICAL, low_memory=False)
    moves = ([(s, v, COL[a], f, t) for s, v, a, f, t, *_ in UPDATES]
             + [(s, v, COL[a], "unknown", t) for s, v, a, t, *_ in ADDITIONS]
             + [(s, v, COL[a], f, "unknown") for s, v, a, f, _w in REMOVALS])
    for system, variant, col, frm, to in moves:
        m = (base["system"] == system) & (base["variant"] == variant)
        if m.sum() != 1:
            raise SystemExit("ABORT: %s/%s matched %d canonical rows"
                             % (system, variant, m.sum()))
        got = _cell(base.loc[m, col].iloc[0])
        if got != frm:
            raise SystemExit("ABORT: canonical %s/%s %s is %r, expected %r"
                             % (system, variant, col, got, frm))
        base.loc[m, col] = to
    print("  canonical: %d expectation cells moved" % len(moves))
    out = rederive(base)

    # rederive() grades every row and does NOT consult
    # mech_graded_excluded_reason, so variants excluded from grading by
    # apply_brct_pooling_v73 (BRCA1 R1699L and R1699Q, whose mechanism is not
    # observable on the supplied structure) come back graded. The exclusion is
    # recorded in the canonical, so re-apply it here rather than editing the
    # shared rederive path, which earlier corrections also import.
    EXC = "mech_graded_excluded_reason"
    if EXC in base.columns:
        excluded = base.loc[
            base[EXC].notna() & (base[EXC].astype(str).str.strip() != ""),
            ["system", "variant"]]
        keys = set(map(tuple, excluded.values))
        gcols = [c for c in out.columns if c.startswith("mech_consistency_")
                 or c.startswith("nbhd_mech_consistency_")]
        m = out.apply(lambda r: (r["system"], r["variant"]) in keys, axis=1)
        for c in gcols:
            out.loc[m, c] = None
        print("  grading exclusion re-applied to %d variant(s): %s"
              % (int(m.sum()), ", ".join(sorted(out.loc[m, "variant"]))))
    if not dry_run:
        ts = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
        shutil.copy2(CANONICAL, SNAPDIR / ("scored_61var_canonical.pre_v80_%s.csv" % ts))
        out.to_csv(CANONICAL, index=False)
        print("  wrote    -> %s" % CANONICAL.name)
    return out


def check() -> int:
    base = pd.read_csv(CANONICAL, low_memory=False)
    led = pd.read_csv(LEDGER)
    bad = []
    for system, variant, axis, _f, to, etype, _d, _b, cite in UPDATES:
        r = base[(base.system == system) & (base.variant == variant)]
        if _cell(r[COL[axis]].iloc[0]) != to:
            bad.append("canonical %s/%s %s != %s" % (system, variant, axis, to))
        lr = led[_match(led, system, variant, axis)]
        if len(lr) != 1 or str(lr.expected_token.iloc[0]) != to:
            bad.append("ledger %s/%s/%s != %s" % (system, variant, axis, to))
        elif str(lr.evidence_type.iloc[0]) != etype or cite not in str(lr.evidence_citation.iloc[0]):
            bad.append("ledger %s/%s/%s basis metadata not updated"
                       % (system, variant, axis))
    for system, variant, axis, to, *_ in ADDITIONS:
        r = base[(base.system == system) & (base.variant == variant)]
        if _cell(r[COL[axis]].iloc[0]) != to:
            bad.append("canonical %s/%s %s != %s" % (system, variant, axis, to))
        if _match(led, system, variant, axis).sum() != 1:
            bad.append("ledger row missing for %s/%s/%s" % (system, variant, axis))
    for system, variant, axis, _f, _w in REMOVALS:
        r = base[(base.system == system) & (base.variant == variant)]
        if _cell(r[COL[axis]].iloc[0]) != "unknown":
            bad.append("canonical %s/%s %s should be unknown" % (system, variant, axis))
        if _match(led, system, variant, axis).sum() != 0:
            bad.append("ledger row for %s/%s/%s should be gone" % (system, variant, axis))
    for system, variant, axis, keep, _w in REJECTED:
        lr = led[_match(led, system, variant, axis)]
        if len(lr) != 1 or str(lr.expected_token.iloc[0]) != keep:
            bad.append("REJECTED change %s/%s/%s was applied -- it must stay %r"
                       % (system, variant, axis, keep))
    if bad:
        print("FAIL: v8.0 corrections not consistently applied")
        for b in bad:
            print("   " + b)
        return 1
    print("PASS: v8.0 -- %d updated, %d added, %d removed, %d rejected and "
          "still held; ledger %d rows"
          % (len(UPDATES), len(ADDITIONS), len(REMOVALS), len(REJECTED), len(led)))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        return check()
    patch_canonical(a.dry_run)
    print("  next: python scripts/curate_evidence_types.py  "
          "(regenerates the ledger and its summary from the canonical)")
    if a.dry_run:
        print("  dry run -- nothing written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
