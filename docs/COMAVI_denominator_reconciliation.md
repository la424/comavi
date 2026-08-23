# Denominator reconciliation — adjudicating the external review


<!-- COMAVI_CURRENT_RESULTS_BANNER -->
> **Current denominator reconciliation.** The primary population contains 57 mechanism-gradeable variants and 132 applicable structural outputs; the all-row canonical audit contains 133 structural outputs. The current totals are 99/132 for the primary continuity convention and 99/133 when all 61 resource rows are retained. See [`COMAVI_current_results.md`](COMAVI_current_results.md). The detailed adjudication below is retained as a historical pre-A636P record.
<!-- /COMAVI_CURRENT_RESULTS_BANNER -->

Every number below was regenerated from `reference_outputs/scored_61var_canonical.csv`
through the shipped pipeline (`scripts/apply_concordance_v5.py`), not transcribed.
Reproduce with `scripts/verify_denominators.py`.

## Historical adjudication record (pre-A636P correction)

| # | Claim from external review | Verdict | Evidence |
|---|---|---|---|
| 1 | "Every denominator in this paper coincides at 57" is wrong | **UPHELD** | Denominators are 57 variants / 130 axes / 49 tier-carrying / 63 comparisons — see table below |
| 2 | "63 of the variants carry a directly measured free energy" is wrong | **UPHELD** | 63 = *comparisons*, only 15 in-benchmark; 61-row benchmark cannot contain 63 |
| 3 | "44 variants with a measured destabilization" is wrong | **UPHELD** | 44 = comparisons; 13 in-benchmark, 31 external SKEMPI |
| 4 | BRCT rows were appended under a binary correct/incorrect convention | **FALSE** | Shipped weighted rubric re-run row-by-row reproduces all 12 stored grades |
| 5 | BRCT rows lack the uncertainty-gating fields | **FALSE** | All four monomer gate fields populated 12/12; `ddg_monomer_confident` True for all |
| 6 | 92/120 corrects an erroneous aggregate | **MISLEADING** | 92/120 is the interaction-only subset of the correct 99/131 |
| 7 | The primary mechanism denominator is 47, not 57 | **WITHDRAWN BY ITS AUTHOR** | Reversed later in the same transcript |

## The denominators, as they actually are

The paper's claim that all denominators coincide at 57 is false. They differ
because the metrics count different things, which is legitimate — but it must be
stated, not denied.

| Metric | Denominator | Unit |
|---|---:|---|
| Mechanism-consistency | 57 | graded variants |
| Structural agreement (graded population) | 130 | testable axes |
| Structural agreement (all rows retained) | 131 | testable axes |
| Tier calibration / tier × energy | 49 | tier-carrying variants |
| Monomer-fold axis | 27 | axes |
| Complex-fold axis | 25 | axes |
| Binding axis | 32 | axes |
| Tier axis | 47 | axes |
| Predicted–measured comparisons | 63 | comparisons |
| — of which in-benchmark | 15 | comparisons |
| — measured destabilizers ≥ 1 kcal/mol | 44 | comparisons |

## Measured-energy accounting

From `reference_outputs/COMAVI_delta_calibration_points.csv` (63 rows):

| System | Comparisons | Axis | In benchmark |
|---|---:|---|---|
| BRCA1 BRCT (GdmCl unfolding) | 10 | monomer fold | yes |
| Hb tetramer (assembly energetics) | 5 | binding | yes |
| Barnase–barstar (SPR/ITC, SKEMPI) | 24 | binding | no |
| TEM1–BLIP (SPR/ITC, SKEMPI) | 24 | binding | no |

Of the 44 comparisons with |ΔΔG_measured| ≥ 1.0 kcal/mol, **13 are in-benchmark**
(8 BRCT + 5 Hb) and **31 are external**. Every sentence quoting 44 or 63 must say
"comparisons", never "variants in this study".

## Why the BRCT grades carry no 0.5

The review inferred a different grading convention from the absence of partial
grades. Re-running `grade_mechanism_consistency` on all 12 BRCT rows reproduces
every stored grade (7 consistent, 3 inconsistent, 2 N/A). The absence of partials
is a structural property of a single-axis variant, not a convention:

- Partial (Step 6) requires a correctly-called positive axis *plus* a
  false-positive axis. A one-axis variant has no second axis to false-fire.
- Partial (Step 7 tail) requires a missed axis *without* a "No structural effect"
  call. When the only annotated axis is missed, the call is necessarily
  "No structural effect detected", so Step 5 fires first and returns inconsistent.

Step 5 is applied identically to interaction rows. A single-axis variant can only
reach consistent or inconsistent — it either fires its one axis or it does not.

## Per-cohort axis decomposition

Reconciles exactly: interaction + BRCT = all.

| Cohort | tier | monomer | complex-fold | binding | total |
|---|---:|---:|---:|---:|---:|
| All 61 rows | 34/47 | 21/27 | 20/25 | 24/32 | **99/131** |
| Interaction 49 | 34/47 | 14/16 | 20/25 | 24/32 | **92/120** |
| BRCT 12 | 0/0 | 7/11 | 0/0 | 0/0 | **7/11** |

The BRCT rows *are* evaluated by the shipped monomer gate — they contribute 11 of
the 27 monomer axes. Nothing needs backfilling.

## Adopted reporting scheme

Three populations, each with its purpose stated, replacing the false
coincide-at-57 claim:

| Population | n | MC at t = 2.5 | Purpose |
|---|---:|---:|---|
| All gradeable | 57 | 41/57 = **0.7193** | Performance across all applicable structural contexts |
| Interaction subset | 47 | 34/47 = **0.7234** | Performance where all three axes are potentially evaluable |
| BRCT core/fold subset | 10 | 7/10 = **0.7000** | Monomer-fold performance with direct physical validation (ρ = 0.72) |

The 57-variant population is retained as primary. A monomer-only variant has a
mechanism, and grading only its applicable axis is a legitimate mechanism test;
inapplicable partner-dependent axes are omitted from its expected pattern rather
than scored as correct negatives. This is already what the pipeline does.

## Open decision: R1699L and R1699Q

Currently ungraded by curated role. Under the shipped rubric both regrade
**consistent** — expected class `structurally_silent`, pipeline calls no
structural effect. Including them gives 43/59 = 0.7288 (Δ +0.0095).

Recommendation: keep them excluded, and say so explicitly in §2.7 rather than
leaving the exclusion implicit in the cohort file. Both are pathogenic variants
with an intact fold; scoring the pipeline "consistent" for predicting silence
would award credit for exactly the asymmetry the paper concedes is uninformative
about pathogenicity. They belong in the text as fold-intact functional controls.
