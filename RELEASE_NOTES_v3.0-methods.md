# COMAVI v3.0-methods

A major release. Published numbers changed, so a v2.x citation is superseded rather
than refined. Every value below is produced by a committed generator in this
repository and re-derived by the checks that run in CI.

## Why the major bump

The benchmark's ground truth was audited against the primary sources themselves
rather than against its own recorded summaries, and a defect was found in how one
generator selected its input columns. Both moved published results.

| quantity | v2.1-methods | v3.0-methods |
|---|---|---|
| mechanism-pattern agreement at 2.5 kcal/mol | 0.7193 (41/57) | **0.6607 (37.0/56)** |
| graded for whole-variant localization | 57 of 61 | **56 of 61** |
| committed axis expectations | 107 rows | **99 rows** |
| direction-aware axis agreement | 65/85 = 0.765 | **56/77 = 0.727** |
| four-output continuity aggregate | 99/132 = 0.750 | **86/123 = 0.699** |
| prioritization population | 47 (17 mechanism / 30 no lesion) | **46 (21 / 25)** |
| priority-score ROC AUC | 0.943 | **0.834** |
| top-20 recovery | 16/17 | **17/21** |

The raw agreement fell and the margin over the strongest trivial strategy rose, because
that baseline equals the cohort's no-lesion fraction, which shrank. The benchmark became
harder, not worse. The system-cluster 95% interval on the score is 0.573–0.755.

## Ground truth

Five axis expectations were corrected, each verified against the primary article rather
than against a summary of it:

- **KRAS G12D and G12V, binding** — neutral to destabilizing. The prior basis read
  "WT-like RAF1 binding" and cited Hunter 2015, whose Table 2 reports a 4.8- and
  7.3-fold *decrease* in affinity. The basis misstated the source it named.
- **MLH1 H718Y, isolated subunit** — a stability measurement for this exact variant
  exists (Hinrichsen 2013: moderate destabilization, half-life 64% of wild type).
- **SMAD4 I500V, binding** — unbundled from I500T; an interface measurement for this
  variant exists with directional controls (Lindsay 2025).
- **MSH2 G674R, binding** — commitment *withdrawn*. The prior basis attributed a
  co-immunoprecipitation result to a paper that performed none, and the paper that did
  run co-IP studied G674A, a different substitution at the same position.

A sixth proposed change was **refused**. MSH2 C697F would have been flipped to
destabilizing on a source that could not be retrieved, while the existing neutral
commitment is grounded in Ollila 2006, which tested C697F by co-IP and reports the
interaction intact. Overturning a readable source on the strength of an unreadable one
is not a correction. This cost the release its more favourable number.

## Defect fixed: indistinguishability flags read as binding energies

`build_detection_vs_attribution.py` excluded derived columns by the substring
`_disting`, which does not match `_indistinguishable`. All 25 per-partner
indistinguishability flags were therefore admitted as binding energies, and because the
reducer takes the maximum by absolute value, a flag of 1.0 won whenever no real partner
energy exceeded 1.0 kcal/mol.

Consequences, both now corrected: a fabricated binding energy repeated across rows of the
supplementary variant table, and wrong detection and rejection rates at the 1.0 kcal/mol
screening threshold only. The symptom had been visible throughout — two independent
implementations of the same quantity disagreeing by nine variants — with nothing
comparing them. The fix asserts that no returned column is a near-fully-populated
boolean, and the assertion fires on the old rule.

## Supplementary package

Four of five supporting figures shipped under the wrong number after a cyclic
renumbering that the generator absorbed and the package did not; one shipped file was
byte-identical to a different current figure. Every text check passed throughout,
because none of them reads an image. `verify_si_package.py` now gates byte-identity
against current renders, caption agreement across the manuscript, mapping record, index
and README, and advertised file sizes.

## Provenance

Every numeric token in the manuscript and the full supplement was traced to a committed
record. Values that nothing recomputed were given generators rather than being left as
saved tables: the predicted-versus-measured calibration table, the per-axis threshold
test, the sub-threshold evidence test and the two priority-score sensitivity analyses.
A behaviour audit now tests the Methods' claims against the shipped code, which is a
different failure mode from a stale number and was previously ungated. Claims that could
not be tested honestly are listed as untestable with their reason rather than asserted.

Source verification status for all 99 committed axes is published in
`docs/audit/COMAVI_axis_source_verification.csv`, including the axes that rest on an
inaccessible source.

## Supplementary evidence stratification

S1 Text's evidence-stratification paragraph and the record behind it both moved. The
record had not been regenerated since the ground-truth corrections, so the prose traced
cleanly to a stale record and the two agreed with each other while both were wrong — a
pairwise check cannot see that, and the provenance audit did not.

| quantity (S1 Text, Note S2) | v2.1-methods | v3.0-methods |
|---|---|---|
| committed axes / total gradeable | 107 / 131 | **99 / 124** |
| agreement, E1–E2 quantitative evidence | 21/33 = 0.636 | **20/34 = 0.588** |
| agreement, E3–E5 qualitative or inferred | 44/51 = 0.863 | **36/44 = 0.818** |
| crude contrast, Fisher p | 0.0307 | **0.0412** |
| Mantel-Haenszel OR across expected-state strata | 0.711 (p = 0.869) | **0.758 (p = 0.980)** |
| E1 quantitative energetic axes | 18/25 = 0.720 | **18/26 = 0.692** |

The conclusion is unchanged: stratified by expected token the association is absent, so
evidence grade still does not explain the headline agreement. Only the numbers moved.

## Verification infrastructure

Remote CI had been failing on every push for sixteen days and the second workflow had not
run at all for eight weeks, because it triggers only on the default branch and that branch
had not moved. Local verification had been reported as passing throughout: both reports
were accurate about what they measured, and the failing check was not in the hand-picked
list the local runs executed. `scripts/run_ci_locally.py` now extracts every `run:` block
from every workflow and executes it, so the local and remote sets cannot drift apart.

What that surfaced, none of which changed a headline number:

- The canonical table carried 96 literal `N/A` strings written by a regrade. `N/A` is in
  pandas' default missing-value list and nothing here overrides the read or write, so the
  sentinel could not survive a round-trip — which is exactly the contract the shipped
  migration test asserts. Removing them is byte-provably the only change to the table.
- `verify_stage6.py` still held the pre-correction expectations. Its own re-derivation half
  was passing; only the frozen literals were stale, and they were re-pinned from
  independent generators rather than from that script's own output.
- Two verification gates imported packages declared only in `requirements-analysis.txt`,
  which this workflow does not install, so each failed from the day it was added.
- The figure gate reported a runner without Arial as a figure defect, conflating "the
  figure requests Arial" with "this machine can draw Arial". The first is checked
  everywhere; the second, and the byte comparison and label-overlap geometry that depend
  on it, only where the typeface is present — and where it is absent the gate says so
  rather than passing quietly. Figure freshness is therefore not verified on CI; it is
  verified where the figures are produced, and the committed files are still checked
  against the manuscript and the submission package by font-independent gates.
- The release checks filtered on the legacy figure directory, so a change to a submitted
  figure or to its generator did not trigger them at all.
- A historical correction script asserted a value a later, better-evidenced correction had
  overturned. Its check now declares the superseding value and verifies each row against
  whichever is currently authoritative, so it still detects drift on the row that stands.

## Compatibility

- Ground-truth tokens, cohort membership and all derived records changed. Analyses
  pinned to v2.x outputs will not reproduce against this release.
- The scoring pipeline itself is unchanged in its energetics; FoldX is still supplied by
  the user and never bundled.
- `submission/COMAVI_v7.15_SI/` was renamed `submission/COMAVI_v9.2_SI/` to match its
  contents. The gate resolves this directory by glob, so the name is not load-bearing.
