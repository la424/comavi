# COMAVI ISDS-v1 figures

This directory contains the current publication-oriented COMAVI figure suite.

1. `figure1_unified_comavi_workflow.png` — unified COMAVI workflow.
2. `figure2_population_map.png` — benchmark and analysis populations.
3. `figure3_mechanism_localization.png` — mechanism-localization performance.
4. `figure4_isds_definition_performance.png` — ISDS-v1 definition and ranking performance.
5. `figure5_context_components_and_states.png` — structural-context components and evidence states.
6. `figure6_threshold_tradeoff.png` — recovery-versus-rejection operating-point tradeoff.
7. `figure7_alphamissense_isds_mechanism.png` — AlphaMissense, ISDS-v1, and mechanism information.

`ALT_TEXT.json` contains the release alt text. `SHA256SUMS.txt` binds the
figure files, plotting source, and documentation.

Regeneration command:

    python figures/isds_v1/make_figures.py \
      --analysis-dir reference_outputs/isds_v1 \
      --out-dir PATH_TO_OUTPUT_DIRECTORY

The numerical source of truth is `docs/COMAVI_current_results.md`.
