#!/usr/bin/env python3
"""The benchmark population: the 56 graded variants of the 61 scored.

WHY THIS EXISTS
---------------
COMAVI scored 61 candidate variants, and earlier versions reported all 61 as
"the benchmark" while grading only 56. Five candidates have no gradeable
expectation, and keeping them in the headline cohort forced every downstream
count to be explained relative to a population no metric used: 61 -> 56 graded
-> 46 prioritized, with a 47-variant pathogenicity set that kept three of the
ungraded variants and so nested in neither direction.

From v11.0 the benchmark is the 56 graded variants. The five are still scored,
still shipped in scored_61var_canonical.csv (the table is what the pipeline
produced, and its name is part of the released v3.0.1 schema), and are listed
with their reasons in S2 Table, excluded_variants sheet. Every generator that
needs the benchmark takes it from benchmark() below, so membership is defined
once.

The exclusion list is explicit rather than inferred from grading, and
benchmark_mask() asserts that the two agree. A curation change that made one of
the five gradeable, or made a benchmark variant ungradeable, fails here instead
of silently resizing the benchmark.
"""
from __future__ import annotations

import pandas as pd

GRADES = ("consistent", "partial", "inconsistent")

# (system, variant, reason code, one-line reason). Reason text is the curated
# basis, restated; see S2 Table, excluded_variants sheet.
EXCLUDED_CANDIDATES = (
    ("troponin_ic", "R162W", "no_committed_expectation",
     "Increases troponin C affinity, an interface-stabilizing effect reported "
     "without a magnitude the benchmark could commit on a modeled axis "
     "(Elliott et al. 2000)."),
    ("vhl_elonginc", "W117R", "no_committed_expectation",
     "W117 lies in the beta-domain HIF-binding region rather than at the "
     "modeled elongin C interface; no modeled axis is committed "
     "(Stebbins et al. 1999)."),
    ("msh2_msh6", "G674R", "no_committed_expectation",
     "The basis once committed for this variant reported preserved MSH6 binding "
     "by co-immunoprecipitation, but its cited study performs none; the only "
     "co-immunoprecipitation at this position was of MSH2 G674A, a different "
     "substitution, and no study of G674R tests a modeled axis. Commitment "
     "withdrawn in v8.0 (scripts/apply_ledger_corrections_v80.py)."),
    ("brca1_brct", "R1699L", "mechanism_not_observable_on_structure",
     "Disrupts phosphopeptide recognition; the phosphopeptide partner is absent "
     "from the supplied single-chain structure (Rowling et al. 2010)."),
    ("brca1_brct", "R1699Q", "mechanism_not_observable_on_structure",
     "Disrupts phosphopeptide recognition; the phosphopeptide partner is absent "
     "from the supplied single-chain structure (Rowling et al. 2010)."),
)

N_SCORED = 61
N_BENCHMARK = 56


def excluded_keys():
    return {(s, v) for s, v, _, _ in EXCLUDED_CANDIDATES}


def benchmark_mask(canon: pd.DataFrame) -> pd.Series:
    """Boolean mask of benchmark rows, asserted equal to the graded rows."""
    keys = excluded_keys()
    member = pd.Series([(s, v) not in keys for s, v in zip(canon["system"], canon["variant"])],
                       index=canon.index)
    if "mech_consistency_t25" in canon.columns:
        graded = canon["mech_consistency_t25"].isin(GRADES)
        disagree = canon.loc[member != graded, ["system", "variant"]].values.tolist()
        if disagree:
            raise AssertionError(
                "benchmark membership disagrees with grading for %s; update "
                "EXCLUDED_CANDIDATES deliberately, not by accident" % disagree)
    found = {(s, v) for s, v in zip(canon["system"], canon["variant"])} & keys
    if len(canon) == N_SCORED and found != keys:
        raise AssertionError("excluded candidates missing from the scored table: %s"
                             % sorted(keys - found))
    return member


def benchmark(canon: pd.DataFrame) -> pd.DataFrame:
    """The benchmark rows of a scored table (a copy)."""
    out = canon[benchmark_mask(canon)].copy()
    if len(canon) == N_SCORED and len(out) != N_BENCHMARK:
        raise AssertionError("benchmark has %d rows, expected %d" % (len(out), N_BENCHMARK))
    return out


def excluded(canon: pd.DataFrame) -> pd.DataFrame:
    """The scored-but-excluded candidates, with reason code and text."""
    info = {(s, v): (c, t) for s, v, c, t in EXCLUDED_CANDIDATES}
    out = canon[~benchmark_mask(canon)].copy()
    out["exclusion_reason_code"] = [info[(s, v)][0] for s, v in zip(out.system, out.variant)]
    out["exclusion_reason"] = [info[(s, v)][1] for s, v in zip(out.system, out.variant)]
    return out
