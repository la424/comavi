# Supporting-information inputs

`ledger_citation_keys.csv` maps each distinct `evidence_citation` string in
`reference_outputs/COMAVI_evidence_ledger.csv` to the keys of
`docs/manuscript/references.json` that it cites. `scripts/build_supplement.py`
renders these keys with the manuscript's own first-appearance numbering, so the
reference column of S1 Table follows the manuscript when citations move.

The mapping was derived once (v11), not typed: the reference numbers in the
v10.1 S1 Table were decoded with the v10.1 manuscript numbering (commit e1dcbb3,
the revision whose numbering reproduces all 30 citations), and every key was
confirmed against an independent surname-and-year match to references.json.
A citation that appears in the ledger without a row here stops the build.

`whole_variant_negatives.csv` holds the curated basis and citation of the four
non-disruptive variants graded only at the whole-variant level (no axis-level
expectation is committed for them, so the evidence ledger has no row to carry
the text). The rows are the v9.2 S1 Table text without SMAD4 I500V, which gained a
committed binding expectation in v8.0 and is now in the ledger. One phrase was
made exact: CALM1 D96V read "Site pLDDT 69.5 in the model", which named neither
model and rounded 69.45; it now gives both values, which the builder asserts.
`scripts/build_supplement.py` stops if this set differs from the benchmark
variants that have no ledger row.
