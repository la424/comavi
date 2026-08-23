# Post-merge runner output audit

The public repository contains separate benchmark/evaluation, CHD, and live/external entry points. After applying the shared-engine patch:

1. Run a benchmark fixture and save its final CSV.
2. Run the CHD workflow on a small fixture and save its final CSV.
3. Run the live/external workflow on a small fixture and save its final CSV.
4. Verify all three files:

```bash
python verification/verify_isds_result_csv.py \
  results/benchmark_fixture.csv \
  results/chd_fixture.csv \
  results/live_fixture.csv
```

5. Search runner/report code for explicit output-column allowlists and add all ISDS-v1 fields if needed.
6. Confirm JSON/report exports include the score, both components, dominant axis/partner, availability, and version.
7. Confirm adding unrelated variants to an input batch does not change existing ISDS values.
8. Run the complete repository CI suite from a clean checkout.

This audit is required because core integration creates the fields, but a downstream writer can still drop them if it exports a fixed legacy column subset.
