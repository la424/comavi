# COMAVI — Complex-Aware Variant Impact Scoring

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21284083.svg)](https://doi.org/10.5281/zenodo.21284083)
[![verify](https://github.com/la424/comavi/actions/workflows/verify.yml/badge.svg)](https://github.com/la424/comavi/actions/workflows/verify.yml)

COMAVI is a FoldX-based structural variant-interpretation pipeline. Most structural
variant-effect tools score a mutation against a single protein in isolation. But
disease genes act in **protein complexes**, and a variant's real structural
consequence often only appears in that multimeric context — at an interface, or in
the fold of a subunit as it sits within its complex. COMAVI's central premise is that
variant disruptiveness should be evaluated **in the appropriate multimer**, and that
doing so resolves the **specific mechanism** of disruption rather than emitting a
single undifferentiated score.

**Try it in your browser (no install):** the guided Colab notebook helps retrieve,
upload, or optionally predict the required structures, runs COMAVI, and shows
per-variant mechanism cards.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/la424/comavi/blob/main/notebooks/COMAVI_colab.ipynb)


<!-- COLAB_PUBLIC_CLOSEOUT_V1 -->
### Public Colab workflow

The browser notebook runs the generic `run.py` path for user-supplied genes and
missense variants. It reports both the ISDS-v1 priority score and the signed
mechanism profile, validates residue identity before FoldX scoring, supports
explicit per-gene monomer and multimer numbering offsets, and produces a
verified downloadable report bundle. A factory-fresh two-variant HBB smoke test
and its limitations are recorded in
[`docs/COLAB_PUBLIC_WORKFLOW_VALIDATION.md`](docs/COLAB_PUBLIC_WORKFLOW_VALIDATION.md).

## Start here

Choose the path that matches what you want to do:

1. **Run your variants in a browser** — open the Colab notebook above. It is the
   recommended first-run path and requires no local Python installation.
2. **Check a local installation** — follow the [copy-and-paste installation and
   dry-run](#installation). The HBB–HBA1 example validates configuration,
   structures, and variant expansion without invoking FoldX.
3. **Configure your own genes or reuse a prepared system** — use the
   [setup wizard](docs/COMAVI_SYSTEM_SETUP_WIZARD.md) in Colab, or follow
   [Run it on your own genes](#run-it-on-your-own-genes) for the command line.
4. **Browse the published results** — open
   [`reference_outputs/COMAVI_results_summary.xlsx`](reference_outputs/COMAVI_results_summary.xlsx)
   for the per-variant calls and benchmark metrics, or see the indexed
   [`figures/`](figures/README.md) and
   [`reference_outputs/scored_61var_canonical.csv`](reference_outputs/scored_61var_canonical.csv).
5. **Reproduce or audit the study** — start with the
   [documentation index](docs/README.md), the
   [canonical benchmark ledger](docs/COMAVI_v7_canonical_benchmark_ledger.md),
   and the [no-FoldX self-test](#quick-self-test-no-foldx-required).

<!-- COMAVI_SETUP_WIZARD_V1 -->
### Plug-and-play system setup

The Colab notebook includes a sequence-aware setup wizard. It ranks chain
assignments, infers uniform numbering offsets, checks every submitted residue,
and presents a traffic-light preflight before FoldX. A reviewed system can be
saved as `COMAVI_system_setup_bundle.zip`; collaborators can upload one or more
bundles and score new variants without editing Python or YAML. Prepared bundles
are revalidated against each new variant batch and their offsets are rebuilt for
the current numbering convention. See the
[system-setup guide](docs/COMAVI_SYSTEM_SETUP_WIZARD.md) for supported cases and
explicit limits.

## What it does

For each missense variant, COMAVI computes a **three-axis ΔΔG** profile against
supplied monomer and complex structures using FoldX. Complexes may be experimental
structures or predicted models; COMAVI decomposes their structural effect by mechanism:

- `ddg_monomer` — destabilization of the isolated subunit's fold
- `ddg_fold_{partner}` — destabilization of that subunit's fold *within the complex*
- `ddg_binding_{partner}` — disruption of the protein–protein interface itself

This decomposition is what monomer-based tools miss: a variant that looks innocuous
on the lone subunit can be a clear interface disruptor in the assembled complex, and
COMAVI tells you which axis is hit. ΔΔG concordance is **pLDDT-gated** (≥70 strict, ≥50
relaxed) to suppress FoldX artifacts at low-confidence positions, and interface calls
are gated by interface-position pLDDT rather than raw contact count. A **four-way
concordance framework** then integrates the structural tier, FoldX ΔΔG, AlphaMissense,
and Franklin/ClinVar annotations.

**Scope (important, but not the headline):** COMAVI predicts *structural disruption and
its mechanism* — not pathogenicity. Structural disruption overlaps with, but is not
identical to, pathogenicity, so COMAVI reports mechanism and evidence strength
**separately** from phenotype, and "no structural effect detected" is never silently
read as "benign." The benchmark below measures how well its structural calls agree
with literature-grounded structural expectations.

## Repository layout

```
run.py                   generic runner — score ANY genes from a YAML config (main entry point)
scripts/                 core engine + drivers
  comavi_v7/              the COMAVI package (config, foldx_runner, mechanism,
                          evaluation, metrics, concordance, pipeline, ...)
  run_live.py            benchmark driver (live FoldX run)
  apply_concordance_v5.py   four-way concordance (external tools)
  build_report.py        spreadsheet report
  comavi_v7_baseline_correct.py   baseline-correction step (see docs/)
  new_system.py          scaffold a YAML systems block for new genes
  run_chd.py             CHD pipeline driver (Paper 2)
  prepare_chd_input.py   builds CHD variant input (Paper 2)
configs/                 YAML systems configs (chd + benchmark worked examples)
examples/                runnable worked example + its frozen expected output
  hemoglobin_dimer/        config, variants, structure-fetch script, reference_output.csv,
                           compare_to_reference.py (see examples/README.md)
figures/                 the five manuscript figures + figures/README.md index
inputs/
  raw/                   benchmark variant inputs:
                           benchmark_variants_v5.csv  (44-variant PPI set; run_live.py)
                           benchmark_variants_v6.csv  (56-variant set = 44 PPI + 12 BRCA1-BRCT)
  intermediate/          cached intermediates for the no-FoldX self-test
  chd/                   CHD variant inputs (Paper 2)
reference_outputs/       canonical result files:
                           scored_61var_canonical.csv  (61-variant canonical benchmark =
                             v6's 56 + VWF A1-GPIbα and CFH FH1-4-C3b fold-neutral
                             interface disruptors; see ledger §12)
                           COMAVI_results_summary.xlsx  (11-sheet overview) + .docx narrative
                           concordance CSVs, collapsed CHD outputs
data/                    reference inputs (UniProt domain ranges, variant–domain map)
docs/                    user, methods, validation, and maintainer documentation;
                          start with docs/README.md
supplement/              BRCA1-BRCT monomer-fold supplement (measured ΔG data + scoring)
verification/            self-test that reproduces the headline metrics
archive/                 development history — one-off derivation/grading scripts
                          and superseded data; NOT part of the operating pipeline
```

**Not in this repository:** the manuscript sources and the lab-meeting slides are
withheld until preprint/submission. Nothing needed to reproduce the numbers depends
on them — every value they report is regenerated from `reference_outputs/` by the
scripts in `scripts/` and checked by `verification/`. Two audit scripts
(`audit_evidence_claims.py`, `audit_tier_energy_claims.py`) exist to verify that the
manuscript's literals match generated data; without the manuscript present they
report `[SKIP]` and exit cleanly. The claims themselves are documented in
`docs/COMAVI_v7_canonical_benchmark_ledger.md`.

## Installation

COMAVI is developed and verified on Python 3.11. This local smoke test installs
the dependencies, fetches the worked-example structures, and validates the full
configuration path without calling FoldX:

```bash
git clone https://github.com/la424/comavi.git
cd comavi

python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

python examples/hemoglobin_dimer/fetch_structures.py \
  --out examples/hemoglobin_dimer/structures
python run.py \
  --config examples/hemoglobin_dimer/systems.yaml \
  --variants examples/hemoglobin_dimer/variants.csv \
  --structures examples/hemoglobin_dimer/structures \
  --out results/hbb_dry_run \
  --dry-run
```

On Windows PowerShell, activate the environment with
`.venv\Scripts\Activate.ps1` instead of `source .venv/bin/activate`.

The dry run checks dependency installation, configuration parsing, structure
resolution, residue identity, and variant expansion. Full energetic scoring also
requires the following external inputs:

1. **FoldX 5.x** — proprietary, free for academics from https://foldxsuite.crg.eu/.
   Download your own copy and point COMAVI at it:
   ```bash
   export FOLDX_BINARY=/path/to/foldx
   ```
2. **Protein structures** — the command-line runner consumes supplied monomer and
   complex structures. The Colab workflow can retrieve monomers and experimental
   complexes, accept AlphaFold Server uploads, or use optional ColabFold.

For the complete worked example, including full FoldX scoring and comparison with
the frozen reference result, see [`examples/README.md`](examples/README.md).

## The three pipelines

These are the study-reproduction routes. New users scoring their own variants
should normally use the Colab notebook or the
[generic `run.py` route](#run-it-on-your-own-genes) instead. All three study routes
share one engine; they differ only in input and in whether external tools are
layered on.

**1. Benchmark** (61-variant canonical set: 49 PPI across 13 complexes + 12 BRCA1-BRCT)
```bash
export FOLDX_BINARY=/path/to/foldx
python scripts/run_live.py          # -> results/comavi_v7_results.csv (PPI systems)
# BRCT supplement + VWF/CFH fold-neutral interface additions: see ledger §10, §12
# then concordance + evaluation (see scripts/apply_concordance_v5.py --help and docs/)
```

**2. CHD — full (with external concordance)**
```bash
export FOLDX_BINARY=/path/to/foldx
python scripts/run_chd.py                            # structural results
python scripts/apply_concordance_v5.py --help       # then fold in AlphaMissense + Franklin
```

**3. CHD — structural only (multimer structural results, no external tools)**

This is simply the structural stage on its own — run `scripts/run_chd.py` and stop. The
output `results/chd_rerun/chd_structural_results.csv` contains the three ΔΔG axes,
the structural tier, and pLDDT gating, with **no** AlphaMissense / Franklin / ClinVar
columns. Do not run `apply_concordance_v5.py`.
```bash
export FOLDX_BINARY=/path/to/foldx
python scripts/run_chd.py
```

## Run it on your own genes

COMAVI isn't limited to the genes above — it's config-driven. To score variants in your own
proteins, describe your complexes in a YAML file and supply structures; no Python
editing required.

**Check the PDB before you predict.** Where an experimental multimer of your hub and
partner exists it is the better input: the interface is measured rather than predicted,
so no pLDDT confidence gating applies to the binding axis. Several benchmark systems
run on experimental structures for exactly this reason (2HHB, 6XI7, 1JM7, 1SQ0). The
Colab notebook's *Step 2b-0* searches RCSB for you and verifies that candidate entries
actually contain both genes as polymer entities; the same logic is importable as
`notebooks/comavi_helpers.rcsb_find_and_rank(["HUB", "PARTNER"])`. If nothing is
deposited — common for transcription-factor complexes — predict one with AlphaFold.

```bash
# scaffold a config block (prints the YAML + the structure files you need)
python scripts/new_system.py --hub MYGENE --partner PARTNERA --partner PARTNERB

# then run — auto-expands a simple gene,ref_aa,position,alt_aa CSV across systems
export FOLDX_BINARY=/path/to/foldx
python run.py --config my_systems.yaml --variants my_variants.csv \
              --structures ./structures --out results/my_run --dry-run
```

`configs/chd_systems.yaml` is a worked example; `configs/benchmark_systems.yaml` covers the
harder cases (x-ray/NMR, position offsets, non-standard chains, multi-chain complexes). See
the full [adding-your-own-genes walkthrough](docs/adding_your_own_genes.md).


<!-- ISDS_PUBLIC_WORKFLOW_V1 -->
## Verify and report a completed variant run

New benchmark, CHD, live, and generic runs use the same shared engine and
therefore emit the nine versioned ISDS-v1 fields. For a new CHD-focused batch,
use the public wrapper around the generic runner:

```bash
export FOLDX_BINARY=/path/to/foldx
python scripts/run_chd_variants.py \
  --variants INPUT.csv \
  --structures structures \
  --out results/my_chd_batch
```

The wrapper writes `structural_results.csv`, independently verifies the ISDS-v1
fields, and creates a complete augmented table, prioritized table, JSON audit,
and Markdown report under `results/my_chd_batch/isds_v1_report/`.

For arbitrary genes, run `run.py` with a custom YAML config, then apply the same
verification and report tools:

```bash
python verification/verify_isds_output_surfaces.py --require-available \
  results/my_run/structural_results.csv
python scripts/build_isds_variant_report.py \
  results/my_run/structural_results.csv \
  --out-dir results/my_run/isds_v1_report
```

ISDS-v1 is a unitless structural-disruption prioritization index, not a
pathogenicity probability or a validated binary decision rule. Retain its
components, signed mechanism profile, operating point, and model-scope fields.

## Quick self-test (no FoldX required)

The framework's headline metrics can be reproduced from cached intermediates
without running FoldX, AlphaFold structures, or network access. From the
repository root:

```bash
python verification/verify_stage6.py \
  --intermediate inputs/intermediate/comavi_v7_results_with_nbhd.csv \
  --am inputs/AM_variants_comavi_mechanism_test.xlsx \
  --scripts-dir scripts
```

All three arguments are required. Expected result: `11/11` checks `[ OK ]`.
See **`verification/README.md`** for what each check covers.

### The re-deriving track, and why it is the one that matters

Add `--canonical` to re-derive the scored table from the shipped scorer and compare
it against the committed copy, rather than checking the committed copy against
hardcoded expectations:

```bash
python verification/verify_stage6.py \
  --intermediate inputs/intermediate/comavi_v7_results_with_nbhd.csv \
  --am inputs/AM_variants_comavi_mechanism_test.xlsx \
  --scripts-dir scripts \
  --canonical reference_outputs/scored_61var_canonical.csv
```

Expected result: `41/41` checks `[ OK ]`.

The 11-check track is a smoke test. Because it compares stored columns against
expectations written by hand, it passes even when the shipped code can no longer
regenerate the stored table — which is not hypothetical: a gating bug in the tier
axis (a missing-value truthiness error that admitted single-chain BRCT rows into
tier grading) once persisted here while every check reported success. Only the
re-deriving track can detect a scorer/table divergence, so it runs on every push
via `.github/workflows/verify.yml`, together with `verify_denominators.py`,
`verify_tier_construction.py` and `verify_v30_reconciliation.py`.

Run the re-deriving track after any change to the canonical table, the scorer,
the tier formula, or the concordance logic.

## Benchmark results

COMAVI exposes two complementary outputs:

1. **ISDS-v1** ranks variants by the combined strength of continuous energetic
   and structural-context evidence for a modeled structural disruption.
2. The **mechanism profile** preserves the separate monomer-fold,
   complex-context, and binding predictions needed to interpret the proposed
   lesion and choose a follow-up experiment.

ISDS-v1 is a unitless prioritization index. It is not a pathogenicity
probability, a calibrated probability of structural disruption, or an
externally validated classifier. This release does not establish a binary
ISDS cutoff.

The literature-curated resource contains 61 variants across 14 protein
systems. Fifty-seven variants are gradeable for whole-variant mechanism
localization. The 47 tier-carrying interaction variants form the primary
structural-prioritization population.

### Mechanism localization

At the 2.5 kcal/mol reproducibility reference, whole-variant
mechanism-pattern agreement was 0.72 (0.7193, graded n=57), with a weighted total of 41/57.

Direction-aware agreement across energetic axes in the primary 57-variant
population was 65/85 = 0.765:

- monomer-fold: 21/27;
- complex-context: 20/26;
- binding: 24/32.

The historical four-output `structural_agreement` aggregate is retained as a
continuity and denominator audit rather than as the primary mechanism result.
On the same 57-variant population it is 99/132 = 0.750. With all 61 resource
rows retained, it is 99/133 = 0.744.

The all-row decomposition is tier 34/47, monomer-fold ΔΔG 21/28,
complex-fold ΔΔG 20/26, and binding ΔΔG 24/32; denominators sum to 133.
The one additional all-row axis is the monomer-fold axis of BRCA1 R1699Q,
which is retained in the resource but excluded from whole-variant grading by
its curated role.

### Structural-disruption prioritization

On the 47-variant prioritization population, comprising 17 committed modeled
structural mechanisms and 30 variants curated to lack a committed lesion on
the modeled axes:

- ISDS-v1 ROC AUC 0.943;
- average precision 0.818;
- energetic component ROC AUC 0.857;
- structural-context component ROC AUC 0.890;
- system-cluster ROC AUC 95% interval 0.858–1.000.

The ranking reflects the intended limited-budget use case:

- Top 10: 9 structural-mechanism variants;
- Top 20: 16/17;
- Top 25: 17/17.

These are internal, system-aware benchmark results. They do not establish
external generalization or a validated decision cutoff.

### Structural-context tier

The structural-context tier contains no ΔΔG term. In the 49 tier-carrying
interaction variants, the observed pathogenic fractions were Tier 1 100%,
Tier 2 72%, Tier 3 70%, and Tier 4 43% (Spearman ρ = −0.40,
p = 0.0044).

The tier is reported with its components because interface membership
materially contributes to its performance. ISDS-v1, the tier components, the
signed energetic axes, the thresholded mechanism calls, and the complete
mechanism profile should therefore remain visible together.

## Inputs

**Variant CSV** (one row per variant). The first four fields are required for
structural scoring; `system` is optional:

| column | meaning |
|---|---|
| `gene` | gene symbol (matches the system config) |
| `ref_aa`, `position`, `alt_aa` | reference AA, 1-based residue, alternate AA |
| `system` *(optional)* | target a specific configured PPI system; if omitted, the generic runner expands the variant across matching systems |

The full concordance step additionally uses `AlphaMissense`, `AlphaMissense_pathogenicity`,
and `franklin` columns. System → partner/structure mappings are defined in
`scripts/comavi_v7/build_chd_config.py` (CHD) and `comavi_v7/config.py` (benchmark);
adapt these for your own proteins.

**Structures.** The command-line runner expects monomer and complex coordinate
files under the selected structures directory, named as defined in the system
config. Complexes may be experimental biological assemblies or predicted models.
The Colab workflow can retrieve monomers and candidate experimental complexes,
or accept AlphaFold Server and optional ColabFold models.

## CHD data provenance & required acknowledgment

The CHD variants analyzed here derive from the **Gabriella Miller Kids First
Pediatric Research Program (Kids First)**, supported by the Common Fund of the
Office of the Director, National Institutes of Health. The data were obtained
through the Kids First Data Resource Center and dbGaP, and their use is governed
by a Kids First / dbGaP Data Use Certification. Only distinct-variant-level,
de-identified data are included in this repository - no per-individual genotypes,
allele counts, or sample identifiers.

**dbGaP study accession: `phs001138`** — "Kids First Pediatric Research Program in
Congenital Heart Disease," the Kids First subset of the Pediatric Cardiac Genomics
Consortium (PCGC) study `phs001194`.

> **Verify before submission.** Confirm `phs001138` (and its version suffix) is the
> accession named on your own Data Use Certification, and confirm the acknowledgment
> wording the DUC requires. PCGC publishes dataset-specific acknowledgment templates;
> the Kids First CHD template requires the phs accession number to be inserted
> explicitly, and manuscripts not prepared in collaboration with PCGC investigators
> must carry the PCGC non-endorsement sentence in addition to the Kids First statement.

## Reproducing the full study

Code, the worked example, and the result CSVs live here. The full AlphaFold
structure set (large) is best archived separately on **Zenodo/Figshare** with a
DOI linked from this README, which also keeps a citable record of the exact inputs.
(GitHub's Zenodo integration can additionally mint a DOI for this code on release.)

## Caveats

- **Command line:** bring your own structures and licensed FoldX binary.
  **Colab:** use the guided retrieval, upload, and optional prediction routes;
  FoldX is still supplied by the user and is never bundled.
- Variants in disordered / low-pLDDT regions are structurally **unevaluable**
  regardless of substitution severity — this is reported, not silently dropped.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the supported environment, validation
commands, and reproducible bug-report checklist.

## Citation

See `CITATION.cff`.

## License

MIT — see `LICENSE`.
