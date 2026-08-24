# Reusable CHD and external-variant workflow with ISDS-v1

## Which entry point to use

- `scripts/run_chd.py` reproduces the frozen CHD paper cohort.
- `scripts/run_chd_variants.py` evaluates a new variant list against the shipped
  CHD systems configuration and automatically verifies and reports ISDS-v1.
- `run.py` is the generic entry point for new genes or new protein systems.

All three structural paths delegate to the same `comavi_v7.pipeline` engine.
The wrapper does not maintain a second formula.

## New variants in configured CHD genes

Prepare a CSV with `gene,ref_aa,position,alt_aa`. Extra non-identifying columns
are preserved. Then run:

```bash
export FOLDX_BINARY=/path/to/foldx
python scripts/run_chd_variants.py \
  --variants my_chd_variants.csv \
  --structures structures \
  --out results/my_chd_batch
```

The wrapper creates:

- `results/my_chd_batch/structural_results.csv`: complete shared-engine output;
- `results/my_chd_batch/isds_v1_report/chd_with_isds_v1.csv`: complete verified
  table with the nine versioned ISDS fields;
- `chd_prioritized.csv`: compact stable ranking with mechanism fields;
- `chd_isds_summary.json`: provenance and stored-versus-recomputed audit;
- `chd_isds_report.md`: human-readable top-ranked review table.

Use `--dry-run` before FoldX to validate configs, structures, fan-out, residue
identity, and structure provenance. Use `--no-fanout` only when every input row
already has an explicit `system` value.

## Variants in new genes or systems

Create or extend a YAML systems configuration, then call the generic runner:

```bash
python scripts/new_system.py --hub MYGENE --partner PARTNERA
python run.py \
  --config my_systems.yaml \
  --variants my_variants.csv \
  --structures structures \
  --out results/my_run \
  --dry-run
```

After the real run:

```bash
python verification/verify_isds_output_surfaces.py --require-available \
  results/my_run/structural_results.csv
python scripts/build_isds_variant_report.py \
  results/my_run/structural_results.csv \
  --out-dir results/my_run/isds_v1_report \
  --prefix my_run
```

## Historical pre-ISDS results

The report builder can add ISDS-v1 directly when all required tier and energetic
columns are present. Existing columns are preserved. A partial nine-field ISDS
schema is rejected, and a complete stored schema is independently recomputed.
Use `--replace-existing` only after reviewing a reported mismatch.

The public 384-row historical CHD table has been migrated under
`reference_outputs/chd_isds_v1/`. The 144-row
`reference_outputs/chd_concordance_collapsed.csv` table is a downstream
summary and cannot be used to reconstruct ISDS-v1 because its direct energetic
and structural-context fields are absent.

<!-- CHD_TRACKED_MIGRATION_V1 -->

## Interpretation contract

ISDS-v1 is a fixed, cohort-independent structural-disruption prioritization
index. It is not a pathogenicity probability, a calibrated probability of
structural disruption, or a validated binary rule. Review the energy component,
context component, uncapped energy ratio, dominant signed axis, thresholded
mechanism calls, continuous energies, evaluability, and model scope together.

A silent or unavailable result does not establish benignity. It means that the
supplied model and available energetic axes did not support a modeled static
lesion under the stated conditions.

## Data handling

Do not place patient identifiers, sample identifiers, per-person genotypes, or
controlled-access metadata in public result tables or reports. Use distinct,
de-identified variant-level inputs for shared outputs.

## Verify the CHD structure package

The PDB and CIF files are supplied separately from the GitHub repository.
`reference_outputs/chd_isds_v1/chd_structure_manifest.json` binds the
historical CHD package to the filenames required by
`configs/chd_systems.yaml` without exposing local filesystem paths.

Before scoring variants, verify the package:

    python verification/verify_chd_structure_contract.py       --structure-dir /path/to/chd_structures

This checks all 34 required files—21 PDBs and 13 CIFs—against their recorded
sizes and SHA-256 hashes. The normal COMAVI dry run then checks residue
identity, chain and offset configuration, and declared structure provenance.

<!-- CHD_STRUCTURE_CONTRACT_V1 -->
