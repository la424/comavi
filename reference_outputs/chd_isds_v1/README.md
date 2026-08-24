# CHD ISDS-v1 reference migration

This directory contains a deterministic ISDS-v1 migration of
`reference_outputs/chd_concordance_results_FIXED.csv`.

The source table already contained the structural-context tier and the direct
monomer, complex-context, and binding energetic columns required by the fixed,
cohort-independent ISDS-v1 formula. The migration therefore did not rerun
FoldX and did not alter any pre-existing row or column value.

- Source rows: 384
- Source systems: 10
- ISDS-v1 available rows: 101
- ISDS-v1 unavailable rows: 283
- Formula version: ISDS-v1
- Source SHA-256: `00568515ff0b97264679477880c38ee01238129bb554541f8f7bde88d6738460`
- Inferred partner labels: actin, cdh2, ctnnb1, dvl2, gli3, kpna1, kpna6, mdfi, rock2, shroom3, tcf7l1, zic3

`chd_concordance_with_isds_v1.csv` preserves every source column and appends
the nine versioned ISDS-v1 fields. The prioritized CSV is a structural-review
surface, not a pathogenicity ranking or a validated binary decision rule.
Signed energetic values and mechanism calls must remain visible when reviewing
the ranking.

The 144-row `reference_outputs/chd_concordance_collapsed.csv` file is a
historical downstream summary. It lacks the direct energetic and
structural-context fields needed to reconstruct ISDS-v1 and must not be used
as an ISDS source.

Regeneration command:

    python scripts/build_isds_variant_report.py       reference_outputs/chd_concordance_results_FIXED.csv       --out-dir PATH_TO_OUTPUT_DIRECTORY       --prefix chd_concordance       --top-n 50

Verification command:

    python verification/verify_isds_output_surfaces.py       --require-available       PATH_TO_OUTPUT_DIRECTORY/chd_concordance_with_isds_v1.csv

## External structure-package manifest

`chd_structure_manifest.json` contains the path-free filenames, sizes, and
SHA-256 hashes of the 34 PDB/CIF files used by the historical CHD analysis and
validated against the current `configs/chd_systems.yaml` contract. The
structures themselves are not distributed in this directory.

Verify a supplied structure directory with:

    python verification/verify_chd_structure_contract.py       --structure-dir /path/to/chd_structures

<!-- CHD_STRUCTURE_MANIFEST_V1 -->
