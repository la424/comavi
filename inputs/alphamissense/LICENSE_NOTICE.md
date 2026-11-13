# License of the AlphaMissense data in this directory

The MIT license of this repository does not cover the files in this directory.

`AlphaMissense_benchmark_subset.tsv` holds a subset of the AlphaMissense predictions
(Cheng J et al. 2023, Science 381:eadg7492), restricted to the proteins of the
benchmark. The source is the Zenodo record "Predictions for AlphaMissense",
doi:10.5281/zenodo.8208688, file `AlphaMissense_aa_substitutions.tsv.gz`, which is
distributed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0
International license (CC BY-NC-SA 4.0). The subset is distributed under the same
license: https://creativecommons.org/licenses/by-nc-sa/4.0/

Rows are copied verbatim and selected by UniProt accession; only the leading `#` of
the header line is removed. `scripts/fetch_alphamissense_subset.py` rebuilds the
subset from the source record, and `AlphaMissense_benchmark_subset.provenance.json`
records the source URL, the MD5 checksum of the source file as verified against the
record, the accessions kept and the number of rows.
