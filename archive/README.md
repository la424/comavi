# archive/ — development history

Nothing in this folder is part of the operating COMAVI pipeline. These files are
retained for provenance and reproducibility of the published results — they are the
one-off scripts and superseded inputs used to *derive and grade* the benchmark, not
to *run* COMAVI. A user scoring their own variants never touches anything here.

To run COMAVI, see the repository root `README.md` ("The three pipelines" and
"Run it on your own genes"). To reproduce the headline metrics without FoldX, see
`verification/README.md`.

## derivation/ — one-off analysis & grading steps

Post-processing steps applied during development to reconcile and grade the benchmark.
Each was run once against a specific intermediate; their effects are already baked into
the canonical `reference_outputs/` and the v6 CSVs. Ordered roughly by pipeline stage:

| File | What it did |
|---|---|
| `comavi_v7_neighborhood.py` | Pipeline-2 (±3 neighborhood contact) variant. **Tested and rejected** — it degraded the evidence gradient. Kept as the documented negative result. |
| `mech_consistency_plddt_patch.py` | One-time pLDDT reconciliation of the `mech_consistency` metric (raw 0.70 → reconciled 0.73). |
| `relaxed_regrounding_walk.py` | Relaxed-tier re-grounding walk over the contracted axes (graded→unknown reclassification). |
| `tier_ddg_corroboration_patch.py` | Added structural-evidence corroboration columns derived from existing result columns. |
| `apply_chd_concordance.py` | Early CHD-specific four-way concordance assembler. Superseded by the generic `scripts/apply_concordance_v5.py`. |
| `collapse_chd.py` | Surgical per-variant collapse update for the CHD concordance deliverable. |
| `verify_v5.py` | Earlier verification script. Superseded by `verification/verify_stage6.py`. |

## superseded_data/

| File | Note |
|---|---|
| `comavi_v7_results.csv` | Earliest cached FoldX intermediate. Superseded by `inputs/intermediate/comavi_v7_results_corrected.csv` and `..._with_nbhd.csv`, which the self-test actually consumes. |
