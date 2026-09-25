COMAVI supporting information
=============================

Nine supporting-information items. Upload each file below to PLOS Computational
Biology as a separate Supporting Information item, using the item name in the
first column. The captions for these items are listed at the end of the
manuscript file, as PLOS requires; they are not repeated here and must not be
uploaded separately.

Supporting information is published exactly as provided and is not copyedited,
so each file is self-describing: every workbook sheet opens on its own caption
row, and every consolidated document opens on a contents table.

  item       file             what it contains
  ---------------------------------------------------------------------------
  S1 Table   S1_Table.xlsx    6 sheets: benchmark cohort, evidence composition,
                              model scope, context comparators, denominator
                              sensitivity, fold/binding disagreements
  S2 Table   S2_Table.xlsx    9 sheets: operating-point sweep, disagreement
                              classes, tier distribution, evidence states,
                              detection metrics, pathogenicity ranking,
                              robustness, SKEMPI validation, threshold sweep
  S1 Text    S1_Text.docx     2 sections: implementation details and
                              statistical conventions; evidence provenance and
                              denominator sensitivity
  S2 Text    S2_Text.docx     3 sections: sensitivity analyses; mechanism-class
                              decomposition and fold-arm errors; operating-point
                              and calibration diagnostics
  S3 Text    S3_Text.docx     5 sections: pathogenicity comparison; robustness
                              and metric reconciliation; SKEMPI validation;
                              rank-grade overlap; the self-caught ground-truth
                              correction
  S1 Fig     S1_Fig.png       benchmark stress tests (PDF sidecar included)
  S2 Fig     S2_Fig.png       SKEMPI binding-axis validation (PDF sidecar)
  S3 Fig     S3_Fig.png       structural-context comparator performance
  S4 Fig     S4_Fig.png       priority-score transformation and component ranking

Reading this alongside the manuscript you reviewed
--------------------------------------------------

Twenty-nine supplementary items in the reviewed manuscript were consolidated
into these nine. SI_INDEX.csv maps all 29 to their new location: which file,
and which sheet or section within it. Nothing was dropped and no table cell was
altered in the move.

Section numbers inside the consolidated documents are local to their file. The
item name each section carried in the reviewed manuscript is given in that
file's contents table and in SI_INDEX.csv, so a reader holding the earlier
version can locate any passage. Old item designations appear only as
provenance, never as pointers -- every cross-reference in these files resolves
to one of the nine items above.

Figures
-------

The four SI figure files are byte-identical to the output of the generators in
the repository (figures/ and figures/isds_v1/), so any of them can be
regenerated and checked. PDF sidecars are included where the generator emits
one.
