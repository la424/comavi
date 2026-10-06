# COMAVI full provenance audit

What is now mechanically guaranteed, and what is not. Written after a pass that
enumerated every number in the paper and supplement, verified every committed
axis against its primary source, and tested every behavioural claim the Methods
makes against the shipped code.

No published result changed. Every defect below was in a copy, a record nothing
read, a docstring, or a supplementary cell.

## Guaranteed

**Every number in the paper and supplement traces to something.** 2,110 numbers
enumerated across the manuscript and the full supplement, including table cells
and spreadsheet cells. Of the 1,762 that are result values, 1,634 (92.7%) match
a committed record in `reference_outputs/` by exact or rounded value. The
remainder resolve as 113 journal volume and page numbers inside citation
strings, 7 comma-thousands extraction artifacts, 6 values reported by a cited
source, one physical constant, and two derivable from records. Nothing is
unresolved. `COMAVI_number_register.csv` is the row-by-row record.

**No generator emits a hardcoded result.** 764 result-like numeric literals were
swept across 72 scripts and figure generators. Not one executable literal in any
builder or in the shipped figure generator produces a result value: all 62
strong-signature hits in `make_submission_figures.py` are geometry, and every
`n = ` in it interpolates from data. The 66 literals inside verifier assertions
are intentional tripwires.

**The pipeline does what the Methods says, on 14 testable claims.**
`verification/audit_pipeline_behaviour.py` tests each and passes: five BuildModel
replicates per variant, reference values 2.9/2.9/3.5, R as a maximum rather than
a sum, S = (E+C)/2 with C on the tier mapping and E in [0,1), experimental
coordinates exempt from pLDDT gating, the isolated subunit demonstrably not
extracted from the complex (13 variants pair an experimental complex with an
AlphaFold monomer, which extraction cannot produce), 11 AlphaFold plus 3
experimental monomer systems, excluded variants carrying no grade, both
reverse-direction controls under 0.72 kcal/mol (largest 0.714, which is what
sets that bound), the 46 = 21 + 25 prioritization population, mature-chain HBB
numbering, absence of an assay never used as evidence, and exact
Clopper-Pearson intervals.

**45 of 99 committed axes are verified against a value in the primary source,**
26 more against an explicit directional statement. The quantitative-energetic
class is the strongest: 25 of 31 exact. Nothing in any source contradicted a
committed token.

**Four new gates close the classes this audit found,** each with passing
negative controls, all wired into CI. Suite is 38 gates, green from a clean
clone built from the tracked file list alone.

## Not guaranteed

**13 axes rest on sources I could not read** — 11 with no accessible copy, 2
where the source names the variant only by its trivial name (Hb Kansas rather
than N102T). **9 more are population-database citations**, which are a query
rather than an article claim and were not re-queried here.

**6 BRCA1 axes stand at assay level only.** Starita 2015's PMC record is
metadata-only; its abstract confirms the assay measures BARD1 RING binding, but
the per-variant scores are in a supplementary dataset that was not available.

**Five Methods claims have no executable test**, each for a stated reason: FoldX
is licensed separately and cannot be re-run here; the 7-and-7 experimental/
predicted complex split is a manuscript table rather than a derivable column;
the 0.7/0.4 context attenuation is applied upstream and no column preserves the
pre-attenuation score; the per-axis pLDDT admission decisions are summarised
into `site_plddt_status` and not preserved individually; and the permutation
machinery has no independent implementation to compare against.

**Assay-level direction matching by keyword does not work.** An automated scan
of source text near each variant mention raised 13 apparent conflicts with
committed tokens, and all 13 were window artifacts quoting other variants. Only
reading the source settles direction.

**Curated supplementary sheets are not recomputable** — `evidence_ledger`,
`whole_variant_negatives`, `benchmark_variants`, `model_scope`. The derived
sheets are now gated; these are curated content and remain reviewed by hand.

## What this pass actually found

Ordered by how badly each could have embarrassed the paper.

1. **The supplement contradicted itself on a perfect-value claim.** The S3 Table
   `tier_distribution` sheet asserted `0/17 = 0.000` structural mechanisms in
   Tiers 3-4 while the S2 Text prose beside it already said 4/17. Corrected to
   `4/17 = 0.235`, and now gated.
2. **A supplementary workbook cell was unreachable by every gate.** Editing
   `S3_Table.xlsx` to read `11/13 = 0.999` passed `verify_si_package`,
   `audit_evidence_claims` and `verify_denominators`. A checksum manifest cannot
   close this: the stale sheet was already stale at packaging, so its digest
   matched. `verify_si_derived_sheets.py` recomputes the cells instead, and its
   negative control reinjects the original `0/17` defect and catches it.
3. **`discover_partners` returned 377 labels where there are 25 partners,**
   admitting every CI bound and 0/1 flag as though it named a partner. It has 14
   consumers and several use it raw, including the generator behind the headline
   sweep. It corrupted nothing — the full sweep is byte-identical under both
   lists, verified before any change — but it is the same hazard that once did
   produce a fabricated published binding energy.
4. **A record had sat through two ground-truth corrections.**
   `COMAVI_tier_ddg_concordance` was last written at v7.7 and still carried
   pre-correction classes and `n_variants` 47 against the correct 45. Nothing
   read it, no gate touched it, it was absent from CI. Now regenerated with a
   `--check` mode.
5. **Headline values were duplicated in three ungated documents,** stating
   superseded figures as current. Every stale number in this audit was a copy of
   a value that lives correctly elsewhere: the paper is gated, its copies were
   not. The present-tense claims now cite the record; genuine per-revision
   history was reframed rather than deleted.
6. **A stale cohort and a one-sided bound in S2 Text.** The class-mixture
   section carried 25/32 and a fraction of 0.439 with a single upper bound.
   Recomputing gives 29/27, 0.518, and a two-sided interval of 0.418 to 0.588
   that contains the cohort's own fraction, so the conclusion survives and the
   old text understated the constraint.
7. **Version-history narrative in S3 Text**, against the convention that the
   paper presents only its final state. Removed.
8. **One generator docstring stating facts three corrections had falsified.**
   `build_rubric_baselines.py` claimed a score of 0.693 and a cohort of 29/28.
   The values were removed rather than updated and the docstring now points at
   the record, since a prose claim with no generator is what that script exists
   to replace.

## Errors of mine this pass caught

Recorded because they calibrate how much to trust the rest.

- I wrote the class-mixture interval as 0.417-0.589 from my own sweep grid. The
  analysis was already generator-backed and the authoritative values are
  0.418-0.588. Corrected to the record.
- My reverse-control test selected energy columns by prefix and swept up the
  `_ci95_` bounds and flag columns, reporting a largest value of 3.614 — the
  identical column-selection trap that caused defect 3 above, reproduced in a
  brand-new test. Point estimates named explicitly give 0.714.
- Four further behaviour tests failed on first run, all because the test was
  wrong: three names for one axis, `crystal` as the experimental status value,
  and BRCT's null `monomer_structure_type` (it is a monomer-only system whose
  experimental chain *is* its structure).
- My duplicate-headline gate flagged a line already headed "HISTORICAL …
  superseded", because the word "current" appeared two lines above referring to
  something else. An explicit historical label now overrides.
- The step-1 value index treated every file in `reference_outputs/` as
  authoritative, so a value matching only a stale record would have scored as
  traced. Checked afterwards: zero published values matched only the stale
  record, so the hole did not bite — but it is a real limit of the method.
