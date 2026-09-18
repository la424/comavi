# COMAVI ISDS-v1 figures

This directory contains the current publication-oriented COMAVI figure suite.

1. `figure1_unified_comavi_workflow.png` — unified COMAVI workflow.
2. `figure2_population_map.png` — benchmark and analysis populations.
3. `figure3_mechanism_localization.png` — mechanism-localization performance.
4. `figure6_priority_recovery.png` — Fig 6: score distributions and top-k recovery.
5. `figureS4_priority_definition.png` — S4 Fig: energy transformation and component ROC curves.
6. `figureS3_context_comparators.png` — S3 Fig: the three structural-context screens.
7. `figure5_context_components_and_states.png` — repository diagnostic; its evidence-state panel
   duplicates S10 Table and is not a manuscript figure.

`figure4_isds_definition_performance.png` is SUPERSEDED: it was a four-panel composite whose
panels A and C were cited as S4 Fig while B and D were cited as Fig 6. PLOS Computational
Biology requires one file per numbered figure, so it was split into the two files above and
nothing writes it any more. Remove with:  git rm figures/isds_v1/figure4_isds_definition_performance.png
6. `figure6_threshold_tradeoff.png` — recovery-versus-rejection operating-point tradeoff.
7. `figure7_alphamissense_isds_mechanism.png` — AlphaMissense, ISDS-v1, and mechanism information.

`ALT_TEXT.json` contains the release alt text. `SHA256SUMS.txt` binds the
figure files, plotting source, and documentation.

Regeneration command:

    python figures/isds_v1/make_figures.py \
      --analysis-dir reference_outputs/isds_v1 \
      --out-dir PATH_TO_OUTPUT_DIRECTORY

The numerical source of truth is `docs/COMAVI_current_results.md`.
