# ISDS-v1 in CHD and external variant evaluation

ISDS-v1 belongs in the shared `scripts/comavi_v7/pipeline.py` engine, after tier/evaluability calculation and before threshold-specific mechanism classification. This avoids duplicating the formula in benchmark, CHD, live, or external runners.

## Required behavior

Every newly generated result table should preserve these fields:

- `isds_v1`
- `isds_energy_component`
- `isds_context_component`
- `isds_energy_ratio_uncapped`
- `isds_dominant_axis`
- `isds_dominant_partner`
- `isds_dominant_signed_ddg`
- `isds_available`
- `isds_version`

ISDS-v1 is `NA` when no energetic axis or no valid tier is available. It is never imputed as zero.

## CHD migration

Historical CHD tables can be upgraded without rerunning FoldX:

```bash
python scripts/add_isds_v1.py chd_results.csv chd_results_isds_v1.csv
python verification/verify_isds_result_csv.py chd_results_isds_v1.csv
```

## Runner audit after merge

Run at least one benchmark fixture, one CHD fixture, and one live/external fixture. Confirm that each output CSV passes `verify_isds_result_csv.py`, and confirm that any reporting/export code does not select a legacy fixed list of columns that drops the ISDS fields.

The shared-engine patch is sufficient only when the downstream runners write the returned dataframe without dropping newly added fields. If a runner uses an explicit output-column allowlist, add the nine fields above to that allowlist.
