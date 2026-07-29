# verification/ — reproduce the headline metrics without FoldX

This is a supported entry point, not a throwaway. It lets anyone reproduce COMAVI's
headline benchmark numbers from cached intermediates **without** installing FoldX or
re-running the ~24–48 h structural stage.

## Run it

From the repository root:

```bash
python verification/verify_stage6.py \
  --intermediate inputs/intermediate/comavi_v7_results_with_nbhd.csv \
  --am inputs/AM_variants_comavi_mechanism_test.xlsx \
  --scripts-dir scripts
```

All three arguments are **required** — the script has no defaults and will exit with a
usage error if invoked bare. Expected result: `11/11` checks `[ OK ]`.

`verify_stage6.py` consumes the cached intermediate passed via `--intermediate`
(`inputs/intermediate/comavi_v7_results_with_nbhd.csv`),
subprocess-calls the operating concordance + evaluation scripts
(`scripts/apply_concordance_v5.py`, `scripts/run_evaluate.py`), and checks the
reproduced `structural_agreement` and `mech_consistency` at all five thresholds
(t=1.0/1.5/2.0/2.5 + Sapozhnikov) against the stored expected values, printing
`[ OK ]`/`[FAIL]` per check.

(These are the per-threshold refined-framework numbers — e.g. structural_agreement
0.718 and mech_consistency 0.761 at t=2.5 — which roll up to the ledger's headline
`structural_agreement` ≈ 0.77 / mech_consistency 0.70 raw / 0.73 reconciled.)

## Scope

- `verify_stage6.py` — the current, maintained self-test. **Use this.**
- The earlier `verify_v5.py` has been retired to `archive/derivation/` (kept only for
  the audit trail); it is superseded by `verify_stage6.py`.

No FoldX binary, no AlphaFold structures, and no network access are required — the
cached intermediates are committed to the repo.
