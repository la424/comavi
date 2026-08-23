# Integrated Structural-Disruption Score (ISDS-v1)

ISDS-v1 is a fixed, cohort-independent prioritization index that combines COMAVI's strongest normalized energetic signal with its structural-context tier. It does **not** estimate pathogenicity and it does not replace the signed multi-axis mechanism profile.

For each valid energetic axis:

- monomer ratio = `abs(ddg_monomer) / 2.9`
- complex-context ratio = `abs(ddg_fold_partner) / 2.9`
- binding ratio = `abs(ddg_binding_partner) / 3.5`

Let `R` be the maximum valid ratio. The energetic component is:

`E = R / (1 + R)`

For Tier 1, 2, 3, or 4, the context component is:

`C = (4 - tier_number) / 3`

The integrated score is:

`ISDS-v1 = (E + C) / 2`

## Interpretation

- Range: 0 to 1 (approaches 1 but is not a probability).
- Higher values indicate stronger combined energetic and contextual evidence.
- The score is `NA` when no energy axis is evaluable or no valid tier exists.
- No binary ISDS cutoff is validated in the current release.
- The mechanism-call threshold is separate from ISDS; it does not enter the formula.
- Signed raw energies and the mechanism profile remain authoritative for mechanism localization.

## Output fields

- `isds_v1`
- `isds_energy_component`
- `isds_context_component`
- `isds_energy_ratio_uncapped`
- `isds_dominant_axis`
- `isds_dominant_partner`
- `isds_dominant_signed_ddg`
- `isds_available`
- `isds_version`

Because the benchmark, CHD, and live/external pipelines call the shared `comavi_v7.pipeline.run_pipeline` engine, adding ISDS in the core engine propagates the fields to each output. `scripts/add_isds_v1.py` provides a migration path for historical result CSVs.
