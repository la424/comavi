# ISDS-v1 benchmark, CHD, and external-output closeout

## Already complete in the merged engine

The shared pipeline calculates the same nine ISDS-v1 fields for benchmark, CHD,
live, and generic runs. The full-data writers preserve the returned DataFrame,
so newly generated raw CSVs retain the fields.

## Public-output components

The current public-output surface consists of:

1. `scripts/run_chd_variants.py`, the reusable CHD wrapper around the
   configuration-driven `run.py` entry point.
2. `scripts/build_isds_variant_report.py`, which creates the complete augmented
   CSV, prioritized review table, JSON audit, and Markdown report.
3. `verification/verify_isds_output_surfaces.py`, which independently
   recomputes and verifies the nine ISDS-v1 fields.
4. The ISDS-aware display in `scripts/run_evaluate.py`.
5. The benchmark, CHD, and external-variant workflow documented in `README.md`
   and `docs/CHD_ISDS_v1_USER_WORKFLOW.md`.
6. GitHub Actions checks for the synthetic CHD fixture, tracked historical
   migration, public reporting surface, and structure-package contract.
7. `verification/audit_isds_public_surfaces.py`, which checks that benchmark,
   CHD, generic, and live entry points use the shared implementation.

The older `scripts/build_report.py` remains a historical benchmark workbook.
Current benchmark, CHD, and external-variant reports should use the ISDS-aware
report builder.

## Release verification

The public-output implementation is checked by its unit tests, the independent
output verifier, the public-surface audit, the tracked migration regeneration,
and the CHD structure-contract verifier.

## Historical CHD migration

For any stored CHD result table, first try the report builder directly:

```bash
python scripts/build_isds_variant_report.py INPUT.csv \
  --out-dir results/chd_migrated/isds_v1_report \
  --prefix chd
python verification/verify_isds_output_surfaces.py --require-available \
  results/chd_migrated/isds_v1_report/chd_with_isds_v1.csv
```

If a file lacks the required direct energetic or tier columns, it cannot be
reconstructed from a collapsed summary alone; regenerate it from the closest
uncollapsed structural output or rerun the shared engine.

## Completion boundary

The code closeout is complete when the patch is merged and CI passes. A claim
that a historical CHD paper figure or table is synchronized additionally
requires migration or regeneration of its exact source CSV and a source-to-
figure comparison. No patient-level or controlled-access identifiers should be
committed.

## Current-pipeline CHD runtime validation

The runtime closeout used the code and `configs/chd_systems.yaml` from GitHub
main commit `215874ed0721e592fa4160df9c5f2915760546b1`. The follow-up branch did
not modify the shared COMAVI scoring engine.

The 34-file historical CHD structure package matched the current configuration
exactly. A current-runner dry run loaded all 10 systems, expanded three
SHROOM3 variants to 15 system-specific rows, validated 30 monomer/complex
residue positions, and passed structure-provenance validation.

A real one-replicate smoke calculation for `zic3_mdfi / ZIC3 C297F` completed
with zero FoldX failures. It produced all three direct energetic outputs,
an available ISDS-v1 result, the complete structural CSV, independent output
verification, and all four public report artifacts.

This smoke calculation verifies executable integration; it is not a
replacement for the historical multi-replicate result. The one-replicate
energies and thresholded mechanism therefore need not be numerically identical
to the deposited multi-run row.

The historical 203-column result and current 213-column reference share
identical direct energetic outputs. Their schema difference reflects the
legacy `mavis_*` versus current `comavi_*` derived fields and 10 added
path-only concordance fields. A small number of shared external-evidence and
concordance fields also differ; this closeout does not assign a cause to those
external-evidence changes.

The tested macOS FoldX build does not provide a reliable standalone
`--version` response and terminated that probe with an internal assertion.
The real FoldX calculation, rather than the version flag, is the runtime
validation used for this release.

<!-- CHD_RUNTIME_VALIDATION_V1 -->
