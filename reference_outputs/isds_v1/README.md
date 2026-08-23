# ISDS-v1 reference outputs

This directory contains the compact, reproducibly generated outputs for the
fixed, cohort-independent integrated structural-disruption score, ISDS-v1.

The canonical 61-variant benchmark is stored one directory above as
`reference_outputs/scored_61var_canonical.csv`. That table is the public
per-variant source of truth and contains the nine versioned ISDS-v1 fields.

The primary structural-prioritization analysis contains 47 interaction
variants: 17 with committed modeled structural mechanisms and 30 curated to
lack a committed lesion on the modeled COMAVI axes.

ISDS-v1 is a unitless prioritization index. It is not a pathogenicity
probability, a calibrated structural-disruption probability, or an externally
validated classifier. No binary ISDS cutoff is established in this release.

The full 10,000-draw system-cluster bootstrap table is not committed because
the reported intervals and public figures use its compact summary. The full
draw table can be regenerated from the committed inputs.

Regeneration command:

    python scripts/analyze_isds_v1.py \
      --canonical reference_outputs/scored_61var_canonical.csv \
      --tier-calls data/analysis_inputs/isds_v1/tier_comparator_variant_calls.csv \
      --out-dir PATH_TO_OUTPUT_DIRECTORY
