# Contributing to COMAVI

Bug reports, documentation fixes, worked examples, and narrowly scoped code
improvements are welcome. Please open an issue before a large behavioral change
so its scientific and compatibility implications can be discussed first.

## Supported environment

COMAVI is developed and tested with Python 3.11. FoldX is third-party licensed
software and must not be committed to the repository or attached to an issue.

```bash
git clone https://github.com/la424/comavi.git
cd comavi
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Before submitting a pull request

Run the checks relevant to your change from the repository root. The following
commands cover the public contracts and current generated-result summary:

```bash
python -m unittest \
  verification.test_isds_v1 \
  verification.test_isds_output_contract \
  verification.test_update_canonical_isds_v1 \
  verification.test_isds_variant_report \
  verification.test_chd_variant_runner \
  verification.test_chd_migration_reproducibility \
  verification.test_colab_notebook_contract \
  verification.test_colab_setup_wizard

python scripts/build_current_results_summary.py --check
python verification/audit_isds_public_surfaces.py
python verification/verify_colab_notebook_contract.py notebooks/COMAVI_colab.ipynb
python verification/verify_colab_setup_wizard_contract.py
```

Also run `git diff --check`. Pull requests should explain what changed, why it is
needed, which commands were run, and whether any public output or scientific claim
changed. Generated reference outputs should only change when the generating code,
input provenance, and validation record change together.

## Reproducible bug reports

Include:

- the exact COMAVI commit;
- Python and operating-system versions;
- whether the run used Colab or the command line;
- the command or notebook step that failed;
- the complete error message;
- a minimal non-sensitive input example; and
- whether the structure was experimental or predicted.

Do not submit licensed FoldX binaries, access credentials, controlled participant
data, or identifiable clinical information. Redact local paths and tokens from
logs before attaching them.
