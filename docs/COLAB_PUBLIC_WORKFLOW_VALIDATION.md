# COMAVI factory-fresh Colab live-run audit

- Uploaded bundle: `COMAVI_variant_results_bundle.zip`
- Bundle SHA-256: `a89dea7eb64b7a6079cb2ffb01cf6a62257613cc88fb84b3f8f628cfe2325844`
- Files: **12**
- Runtime COMAVI commit: `aeeaa3956b26edd67115083941954727316ca997`
- Formula version: **ISDS-v1**

## Verdict

**PASS.** The factory-fresh Colab run successfully evaluated two distinct user-supplied HBB substitutions through the generic COMAVI runner, after applying the correct HBB monomer numbering offset. The bundle contains complete raw, prioritized, machine-readable, and reader-facing outputs, and the stored ISDS-v1 values recompute exactly.

## Contract checks

- **PASS — zip inventory exact.** 12 files
- **PASS — two input variants.** {('hbb', 'E', 6, 'V'), ('hbb', 'E', 6, 'A')}
- **PASS — two output rows.** {('hbb', 'E', 6, 'V'), ('hbb', 'E', 6, 'A')}
- **PASS — one system.** {'hbb_hba1'}
- **PASS — monomer offset applied.** {-1}
- **PASS — multimer offset applied.** {0}
- **PASS — monomer maps to residue 7.** {7}
- **PASS — complex maps to residue 6.** {6}
- **PASS — residue identity all match.**
- **PASS — structure evaluable.**
- **PASS — multimer evaluable.**
- **PASS — ddg evaluable.**
- **PASS — no FoldX errors.** [nan, nan]
- **PASS — ddg_monomer finite.** [0.5662, 0.1942]
- **PASS — ddg_fold_hba1 finite.** [0.6981, 0.2245]
- **PASS — ddg_binding_hba1 finite.** [0.0041, 0.003]
- **PASS — five monomer FoldX replicates per variant.** [5, 5]
- **PASS — nine ISDS fields present.** set()
- **PASS — ISDS available 2/2.** [True, True]
- **PASS — ISDS version current.** {'ISDS-v1'}
- **PASS — priority order E6V then E6A.** ['E6V', 'E6A']
- **PASS — report recomputation exact.** {'all_fields_match': True, 'failing_fields': [], 'maximum_numeric_delta': 0.0, 'tolerance': 1e-08}
- **PASS — provenance commit.** {'comavi_commit': 'aeeaa3956b26edd67115083941954727316ca997', 'isds_version': 'ISDS-v1', 'input_variants': '2', 'output_rows': '2', 'isds_available_rows': '2'}
- **PASS — provenance row counts.** {'comavi_commit': 'aeeaa3956b26edd67115083941954727316ca997', 'isds_version': 'ISDS-v1', 'input_variants': '2', 'output_rows': '2', 'isds_available_rows': '2'}
- **PASS — run log contains Residue-identity check passed: 4 axis position(s).**
- **PASS — run log contains Structure-provenance check passed.**
- **PASS — run log contains Total FoldX failures: 0.**
- **PASS — run log contains ISDS-v1 available: 2/2.**
- **PASS — run log contains Done: 2 rows.**
- **PASS — stderr empty.** 0 bytes
- **PASS — augmented table preserves all raw columns.** []
- **PASS — augmented table preserves raw values.** numeric=[], text=[], max_delta=0.0

## Result summary

- **HBB E6V:** ISDS-v1 `0.26367620`, dominant modeled axis `complex_context` (`0.6981` kcal/mol), `Tier 3`, mechanism call at 1.0 and 2.5 kcal/mol: `No structural effect detected` / `No structural effect detected`.
- **HBB E6A:** ISDS-v1 `0.03592575`, dominant modeled axis `complex_context` (`0.2245` kcal/mol), `Tier 4`, mechanism call at 1.0 and 2.5 kcal/mol: `No structural effect detected` / `No structural effect detected`.

The HBB E6V result is a software smoke test, not a positive-mechanism control. The supplied 2HHB tetramer does not represent the intermolecular sickle-hemoglobin polymer fiber interface, and residue E6 is not an alpha-beta interface residue in this configuration. Therefore, the absence of a thresholded COMAVI lesion here must not be interpreted as evidence against E6V pathogenicity.

## Notebook improvement identified

The initial run correctly stopped on a residue-numbering mismatch. The successful rerun required `monomer_offset: -1` for HBB because the AlphaFold/UniProt monomer includes the initiating methionine, whereas the historical beta-globin variant and 2HHB use mature-chain numbering. The public notebook should expose per-gene monomer and multimer offset overrides as normal form fields rather than requiring a custom cell.

## File hashes

- `COMAVI_run_provenance.txt` — `3a801bf64abe2974945a637227397de6e4befc40607ea9d7e6b710b26c4e01eb`
- `comavi_mechanism_summary.csv` — `ad6698048139bef884dc945867576e7f3088fc8994fc2ce1afebe78268efe681`
- `comavi_variants_isds_report.md` — `73e462ab4fc6019b77e1b03dfd4d37fb509bee63da5b7d83ab0ca2f4bc3b6b4d`
- `comavi_variants_isds_summary.json` — `14c7eeb42a86e53321a653a1566da41e719b08b14acf0d0f27d41d968dffff7b`
- `comavi_variants_prioritized.csv` — `8fbfef9d32093a8b2ada4dbcea4904d94eafd86d99636977780943ec4a8c3e00`
- `comavi_variants_with_isds_v1.csv` — `d4f0906138a5dab36d9e79e71224b054240afea2893aefb409ac6a464f074f8f`
- `config.yaml` — `f1ca6fdbc86b06f2c11ce2679dbf702376d03fbfacd975ff6ddadbd18ed962b2`
- `per_variant_axes.png` — `0d3644de8e6f763517d46a172560ec8bbe3cf721b500397eaa23032b442d703d`
- `run_stderr.log` — `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
- `run_stdout.log` — `44528275a8955a70ee95b81785c0e6781139ef253f22aa93fdbfe39de053ea06`
- `structural_results.csv` — `28a9b15f5bbc4a76164bd0a7d366b2d132869b2836acbd937faa79107f80c748`
- `variants.csv` — `2cdd1d45c992d5193664e07e4491a9a7ddf8fbb19199f08c70b349f3ea104970`
