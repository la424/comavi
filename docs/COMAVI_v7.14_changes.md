# COMAVI v7.14: what changed and why

v7.14 is a full review of the v7.13 submission draft covering accuracy, citations, structure, prose, figures, tables
and PLOS Computational Biology compliance. The scientific claims are unchanged. Two published numbers were computed
wrongly and have been corrected (sections 1 and 2), several disclosures have been added, and the paper is
restructured and shortened.

## 1. Accuracy

Each quantitative claim was recomputed from its committed record (`COMAVI_v714_accuracy_audit.csv`, 74 rows).
51 recomputed exactly. The others are listed below, and every one is fixed.

- **Measured destabilizations recovered** compared absolute values, so it counted eight sign-inverted predictions as
  recoveries and included two measured *stabilizing* effects. Methods require the sign to be right. With the
  generator corrected (`scripts/build_measured_effect_tables.py`) the values are 27, 23, 18, 12 and 11 of 42 across
  the thresholds; they were 35 to 12 of 44.
- **The SKEMPI check did not disclose an exclusion.** Seven barnase Glu73 substitutions were removed before computing
  ρ = 0.55 (n = 41). All 48 now lead (ρ = 0.34, p = 0.019), and the exclusion is reported as a sensitivity view. The
  exclusion was introduced with the analysis, so nothing describes it as prespecified.
- **The pooled calibration had no generator.** `scripts/build_delta_calibration_stats.py` now writes it; slope 0.42
  and ρ = 0.48 over all 63 comparisons, or 0.43 and 0.63 without Glu73. The slope that had been typed into a figure by
  hand is gone.
- **Table 1** listed no pathogenic VWF variant (R1334Q is pathogenic) and showed 0 for the BRCT system, whose variants
  carry no clinical label (now "–").
- Captions that misdescribed their figures were corrected: Fig 4D called three of the five hemoglobin variants
  external, but all five are benchmark variants; the Fig 5 caption named colours the figure does not use; panel B of
  the old AlphaMissense figure plotted the energy component, not the quantity its caption named; and the
  "W37 series" label included N102T.
- Three "prespecified" and "before results" phrases that the history does not support were removed. Only the
  thresholds keep "fixed before scoring", which the first commit supports.
- Two different 47-variant sets were treated as one. The AlphaMissense comparison includes TNNI3 R162W and VHL W117R
  and lacks VWF R1334Q and A1381T. The text now says so.
- Five whole-variant no-lesion variants (CALM1 D96V, N98S, F142L; SMAD4 I500T, I500V) carry no committed axis. They
  stay in every metric and are now disclosed, with the score without them (36.0/52 = 0.692) as a sensitivity view.

## 2. Citations

The reference list grows from 56 to 75 entries, each checked against CrossRef (`COMAVI_v714_reference_audit.csv`).

- Four entries listed authors who are not on the paper.
- A Discussion sentence attributed FoldX's large-effect compression to two papers that do not report it. It now cites
  the two FoldX evaluations that do, both read in full.
- Eleven sources behind the ground truth were missing from the bibliography and are now cited.
- CALM1 N98S was cited to a paper that only mentions it; the primary report (Nyegaard et al. 2012) is now cited.
- The evidence ledger attributes PMID 11735257 to "Kimura 2001". PubMed gives Takahashi-Yanaga et al. 2001, and the
  S3 Table reference column resolves the citation to that entry.

## 3. Structure and prose

- PLOS order (Introduction, Results, Discussion, Methods). Each idea has one home: `COMAVI_v714_claim_map.csv` maps
  68 ideas to their place, and 22 that appeared in two to nine places were consolidated.
- Everything before the References is 7,398 words, including tables and captions; the v7.13 body alone was 10,847.
  A content ledger compared every decimal, fraction, count, variant and item reference before and after, and each
  value that disappeared was either restatement or was reinstated.
- The manuscript is now built from `docs/manuscript/COMAVI_manuscript_source.txt` by `scripts/build_manuscript.py`,
  which numbers citations in order of first appearance and has a `--check` mode, so the .docx cannot drift from its
  source.

## 4. Figures and tables

- Seven main figures become five, plus five SI figures (`reference_outputs/figure_mapping/COMAVI_v714_figure_lineage.csv`).
  All ten come from `figures/submission/make_submission_figures.py`, drawn at printed width in Arial with every label
  at 8 to 12 pt, as 600 dpi LZW TIFFs for upload and 300 dpi PNGs for the manuscript. Fig 1 is the author's PowerPoint
  layout.
- The generator refuses a figure with a title in the image, the wrong font or size, or overlapping labels. It asserts
  every plotted summary against its record, and both builds are byte-reproducible. Two defects were fixed on the way:
  libtiff wrote one uninitialized pad byte, and separate renders gave the PNG and TIFF different crops.
- New gate `scripts/verify_figure_specs.py`, wired into CI, checks format, resolution, size and the checksum
  manifest. It has negative controls.
- Tables are editable Word tables with one style and a 9 pt font, each placed after the paragraph that first cites it,
  with the title above and the legend below. All 19 captions pass one template (`COMAVI_v714_caption_audit.csv`),
  which includes a check that the panel letters in each legend match the panels drawn.

## 5. Supporting information

Eleven items: S1–S3 Table, S1–S3 Text and S1–S5 Fig, in `submission/COMAVI_v7.14_SI/`.

- S3 Table publishes the gated evidence ledger (98 committed expectations) with a reference column that resolves every
  citation to the main list. A new sheet holds the five whole-variant negatives.
- Curation-process notes ("LABEL CORRECTION", "Pipeline should call…") were rewritten as evidence statements, and
  every earlier-draft item name was removed.

## 6. PLOS compliance

`COMAVI_v714_plos_format_audit.csv`: 25 of 28 requirements pass. The other 3 are for the authors: the author list
and affiliations, the acknowledgments, and the funding, competing-interest and data statements entered in the
submission system. Two gaps were fixed in this pass. BRCT is now expanded at first use, and energies stay in
kcal/mol, the field convention, with the kJ/mol conversion stated.

## Open items for the authors

1. Crotti et al. 2013 reports D96V in CALM2. The benchmark labels it CALM1, which encodes the same protein, so no
   result changes.
2. Is the ledger's "Frazier 2008" the *Pediatr Cardiol* paper?
3. Which gnomAD version was used?
4. The name form "Sánchez-Izquierdo Besora P".
5. Should the ledger's "Kimura 2001" citation string be corrected to Takahashi-Yanaga et al. 2001? That changes the
   ledger file itself, so it needs its checksum pin updated.
