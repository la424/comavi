# verification/ — reproduce the headline metrics without FoldX

This is a supported entry point, not a throwaway. It lets anyone reproduce MAVIS's
headline benchmark numbers from cached intermediates **without** installing FoldX or
re-running the ~24–48 h structural stage.

## Run it

```bash
python verify_stage6.py
```

`verify_stage6.py` consumes a cached intermediate in `../inputs/intermediate/`
(`mavis_v7_results_with_nbhd.csv`; pass a different one with `--intermediate`),
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
