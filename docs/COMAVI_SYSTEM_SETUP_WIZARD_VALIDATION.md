# COMAVI system-setup wizard validation record

## Development status

This package is an implementation candidate for the next public Colab release.
It is not yet the merged GitHub version.

## Automated validation completed

The candidate passed:

- Python compilation of the setup module and verification tools;
- notebook JSON and Python-syntax validation;
- static notebook-contract verification;
- eight notebook mutation tests that require current public behavior and reject
  stale, incomplete, or unsafe variants;
- nine setup-module unit tests covering:
  - historical HBB numbering inference;
  - automatic two-gene chain assignment;
  - distinct monomer/multimer offset derivation;
  - wrong-chain rejection;
  - non-uniform-numbering rejection;
  - traffic-light preflight generation;
  - setup-bundle hash verification, reference-sequence merging, and multi-bundle merging;
  - current-batch revalidation of prepared structures with rebuilt offsets;
  - config-derived, system-scoped run provenance for prepared and multi-system runs;
  - plain-language result cards without probability claims;
- an independent synthetic setup-wizard contract test.

## Existing live evidence retained

The previous public notebook was run in a factory-fresh Google Colab runtime
with a licensed Linux FoldX binary, HBB E6V, synthetic HBB E6A, AlphaFold HBB
and HBA1 monomers, and experimental complex 2HHB. That run established:

- successful generic `run.py` execution for two substitutions;
- successful monomer, assembled-complex, and partner-binding calculations;
- all nine ISDS-v1 fields;
- a verified downloadable result bundle;
- correct failure on a residue-numbering mismatch before FoldX;
- successful rerun after distinct monomer and multimer offsets were supplied.

That live run motivated the automatic sequence and numbering wizard. It does
not by itself certify the new automatic mapping code.

## Live validation required before broad release language

Run the candidate in fresh Colab sessions for this matrix:

1. experimental heterodimer with a known numbering offset;
2. experimental heterodimer with variants in both genes;
3. AlphaFold Server heterodimer;
4. optional ColabFold heterodimer;
5. complex with two distinct partners;
6. truncated experimental construct with a uniform offset;
7. deliberately non-uniform numbering that must stop;
8. low-confidence or uncovered variant that must abstain or stop;
9. prepared-bundle reuse with a new variant list;
10. multi-system run from two setup bundles.

Each case should verify chain assignment, residue identity, configuration,
expected energetic columns, all nine ISDS-v1 fields, provenance, and bundle
contents.

## Claims supported now

- The software architecture is generic across configured systems.
- The candidate automates chain assignment and uniform offset inference.
- The candidate fails safely on wrong chains and non-uniform numbering in
  synthetic tests.
- Prepared setup bundles provide a no-code reuse path, recheck the current
  variants rather than blindly trusting stored offsets, and can be merged for
  multi-system batches.

## Claims not yet supported

- universal new-system accuracy;
- fully automatic biological-state selection;
- one-click support for homomers or repeated gene copies;
- usability without any structural judgment for every protein system;
- prospective scientific validation of ISDS-v1 on new systems.

## Runtime dependency

The setup-wizard module uses `gemmi` for PDB/mmCIF parsing and residue-level
mapping. `gemmi` is therefore a tracked COMAVI runtime dependency and must be
installed in local, continuous-integration, and Colab environments before the
wizard contract tests are run.

<!-- SETUP_WIZARD_GEMMI_DEPENDENCY_V1 -->
