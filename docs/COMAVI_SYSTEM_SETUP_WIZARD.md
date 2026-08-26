# COMAVI system-setup wizard

## Purpose

The setup wizard turns the public Colab notebook from a manually configured
structural workflow into a guided, safety-checked system builder. It is meant
to make two common use cases straightforward:

1. **Add new variants to an existing prepared system.** Upload a validated
   `COMAVI_system_setup_bundle.zip`, enter or upload variants, and run COMAVI.
2. **Prepare a new heteromeric system.** Enter a hub gene and partner genes,
   retrieve or upload monomers, review ranked experimental complexes or upload
   a predicted complex, confirm the preflight, and save a reusable bundle.

The wizard automates technical mapping. It does not replace biological judgment
about whether a structure contains the relevant ligand, modification,
conformation, oligomer, or disease state.

## What is automated

### Variant input

The notebook accepts either:

- inline entries such as `SHROOM3 G1003R`; or
- a CSV containing `gene`, `ref_aa`, `position`, and `alt_aa`.

The four required fields are normalized and validated before structural work.
Extra columns are preserved.

### Reference-numbering reconciliation

For each gene, the wizard compares submitted reference amino acids with the
reference sequence. It infers a uniform offset under the COMAVI convention:

```text
reference_position = submitted_position - input_offset
```

A unique solution is green. Multiple equally compatible solutions are yellow
and must be reviewed. No common solution is red.

### Structure-chain mapping

The wizard extracts observed amino-acid sequences from every coordinate chain,
then aligns each reference sequence to every chain. Chain assignment uses:

- sequence identity;
- reference and chain coverage;
- coverage of submitted variant positions;
- numbering uniformity; and
- separation from alternative chain assignments.

The assignment is global: one configured gene is matched to one chain. Manual
chain overrides remain available for expert correction.

### COMAVI numbering offsets

For a mapped chain, the wizard estimates the dominant relationship between
reference positions and coordinate residue numbers. It combines that
relationship with the submitted-numbering offset to derive the integer offset
used by COMAVI.

The current COMAVI engine requires one uniform integer offset per gene and
structural context. Structures with internal renumbering discontinuities are
red rather than silently coerced.

### Experimental-structure ranking

The Colab notebook can download a small RCSB candidate set and assess each
candidate using the mapping criteria above. Candidate ranking also considers
reported resolution when available.

The top candidate is a recommendation for review, not an automatic biological
assertion. The notebook requires an explicit confirmation that the chosen
coordinates represent the relevant partners, construct, ligand/cofactor state,
conformation, modification, and oligomer before a new-system bundle can run.

### Traffic-light preflight

Every submitted variant receives a preflight row:

- **READY / green** — reference amino acid and numbering map in both monomer and
  complex contexts;
- **REVIEW / yellow** — technically usable but an ambiguity or weaker mapping
  needs explicit confirmation;
- **STOP / red** — missing residue, wrong chain, non-uniform numbering, or
  another condition that could mutate the wrong position.

The generic `run.py` entry point repeats residue-identity and structure-source
checks before FoldX. The wizard adds guidance; it does not weaken the engine's
fail-safe behavior.

## Prepared setup bundles

A successful new-system preflight can produce:

```text
COMAVI_system_setup_bundle.zip
```

The bundle contains:

- `config.yaml`;
- required coordinate files;
- `preflight.csv`;
- `setup_report.json`;
- reference sequences;
- a hash and size manifest.

The bundle is verified on upload. For every new batch, the notebook reloads the
reference sequences, rechecks each submitted reference amino acid against the
stored monomer and complex, and rebuilds the combined COMAVI offsets for the
numbering convention used in that batch. A prepared bundle is therefore not a
blind cache of the original offsets.

Users can upload more than one setup bundle; the notebook merges their
configurations, reference sequences, and structures, allowing one variant CSV
to span multiple prepared systems.

Prepared bundles are the recommended route for collaborators with limited
bioinformatics experience. A structurally experienced user prepares and
reviews each system once; downstream users can then add variants without
editing Python or YAML.

## FoldX handling

FoldX remains third-party licensed software and is not distributed by COMAVI.
The notebook supports:

- uploading the Linux x86-64 binary for the current session; or
- selecting a binary stored in Google Drive.

The binary header is checked before use. macOS and Windows builds are rejected
with a direct explanation.

## Result presentation

The notebook retains all raw fields and adds plain-language review cards that
show:

- ISDS-v1 value and rank within the submitted set;
- dominant modeled structural axis;
- threshold-specific mechanism call when available;
- an experiment class suited to that axis; and
- the warning that ISDS-v1 is not a probability or pathogenicity verdict.

The detailed mechanism cards, signed energies, tier, uncertainty, public CSVs,
JSON summary, Markdown report, plots, logs, configuration, preflight, setup
bundle, and exact COMAVI commit remain downloadable.

## Current support boundary

### Supported by the wizard

- new substitutions in a prepared system;
- heteromeric systems with one configured gene per chain;
- variants in the hub or configured partners;
- experimental PDB, AlphaFold Server, and optional ColabFold complex routes;
- distinct monomer and multimer numbering conventions;
- several prepared systems in one batch;
- explicit manual overrides when automatic mapping is reviewable but not
  decisive.

### Requires expert handling or future development

- symmetric homomers or repeated copies of the same gene;
- structures whose numbering cannot be represented by one integer offset;
- isoform choice when gene-symbol resolution is biologically ambiguous;
- deciding whether a coordinate model represents the relevant functional
  state;
- membrane, ligand, post-translational, polymeric, or dynamic mechanisms not
  represented by the supplied static model;
- completely unsupervised selection among biologically distinct assemblies.

## Recommended public wording

Use:

> COMAVI provides a guided browser workflow for new variants and for new
> heteromeric systems. Technical sequence, chain, residue-coverage, and
> numbering checks are automated; biologically important structure choices
> remain explicit and reviewable.

Avoid:

> Upload any variant and COMAVI will automatically determine its biological
> mechanism.

The latter claim is not supported by a static structural pipeline.
