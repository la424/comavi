# COMAVI manuscript — operative writing rules

Distilled from four sources, in priority order. These are the rules the rebuild is
written against. Cited by rule number in review comments.

Sources:
- **[GS]** Gopen GD, Swan JA. *The Science of Scientific Writing.* American Scientist
  78:550–558 (1990).
- **[W]** Whitesides GM. *Whitesides' Group: Writing a Paper.* Adv Mater 16:1375–1377 (2004).
  DOI 10.1002/adma.200400767
- **[TSR-S]** Mensh B, Kording K. *Ten Simple Rules for Structuring Papers.*
  PLOS Comput Biol 13:e1005619 (2017).
- **[TSR-W]** Weinberger, Evans, Allesina. *Ten Simple (Empirical) Rules for Writing
  Science.* PLOS Comput Biol 11:e1004205 (2015). / Mensh & Kording et al., *Ten Simple
  Rules for Writing Research Papers*, PLOS Comput Biol 10:e1003453 (2014).
- **[TSR-F]** Rougier, Droettboom, Bourne. *Ten Simple Rules for Better Figures.*
  PLOS Comput Biol 10:e1003833 (2014).

---

## A. Structure of the whole (outline before prose)

**A1 [W].** The paper is an organized description of hypotheses, data, and conclusions.
Write the *outline* — which includes the data, as finished figures and tables — and settle
it before writing prose. An outline contains little text. Text is cheap; organization is
the expensive part, and it is what must be agreed first.

**A2 [W].** Organize the paper around assimilable data objects — figures, tables,
equations — rather than around text. The text explains the data and is secondary. The more
that is carried by figures and tables, the shorter and more readable the paper.

**A3 [W].** Organize in order of importance, **not chronologically**. Do not recite the
experimental program from initial failures to a climactic finale. Start with the most
important result and put secondary results later, if at all. The reader does not care how
the result was arrived at, only what it is.
→ *For COMAVI:* the benchmark's development history (v5 → v7, the neighborhood pipeline
that was tested and rejected, the 44 → 56 expansion) belongs in supplement or nowhere.

**A4 [W].** Section headings must be specific and information-rich — a claim, not a topic.
"The Rate of Self-Exchange Decreases with the Polarity of the Solvent", not "Measurement of
Rates".
→ *For COMAVI:* "3.4 No benign variant destabilizes on any axis", not "3.4 Axis results".

**A5 [W].** Conclusions are not a summary. They add a higher level of analysis and state
the significance explicitly.

**A6 [TSR-S].** One paper, one central contribution. Everything in the paper exists to
support that single contribution; anything that does not, is cut or moved to supplement.

**A7 [TSR-S].** Context–Content–Conclusion at every scale: the paper, each section, each
paragraph. Open with what the reader needs to orient; deliver the substance; close with
what it means. Never end a unit on a loose end.

**A8 [TSR-S].** The Results are a *sequence of statements*, ordered so each one logically
follows the last, and together they support the central contribution. Order by logic, not
by chronology or by which analysis was run first.

**A9 [TSR-S].** Allocate disproportionate effort to the title, the abstract, and the
figures. Most readers read only these.

---

## B. Structure of the sentence and paragraph (Gopen & Swan)

The seven GS principles, verbatim in substance:

**B1.** Follow a grammatical subject as soon as possible with its verb. Anything of length
between subject and verb is read as an interruption and therefore as less important —
regardless of what it actually contains.

**B2.** Place in the **stress position** (the point of syntactic closure — the end of the
sentence, or the material before a colon/semicolon) the new information you want the reader
to emphasize. Readers exert maximal emphasis on what arrives last.

**B3.** Place the person, thing, or concept whose *story* the sentence tells at the
beginning of the sentence — the **topic position**. Readers expect a unit of discourse to be
a story about whoever shows up first. ("Bees disperse pollen" and "Pollen is dispersed by
bees" are both good sentences; which one is correct depends on whose story the paragraph is
telling. Passive voice is not the defect — a mismatched topic is.)

**B4.** Place appropriate **old information** in the topic position, for linkage backward
and context forward. GS: misplacement of old and new information is the single most common
defect in professional writing.

**B5.** Articulate the action of every clause in its **verb**. If the list of verbs in a
paragraph is `is / are / has / are presumed to be`, the actions are missing and the reader
must guess them.

**B6.** Provide context before asking the reader to consider anything new.

**B7.** Ensure the relative emphasis of the substance coincides with the emphasis the
structure creates.

**B8 [GS].** A sentence is too long **when it has more viable candidates for stress
positions than there are stress positions available** — not at any fixed word count. A
semicolon or colon creates a second stress position and legitimately extends a sentence.

**B9 [GS].** These are principles, not rules. They may be violated deliberately for
effect — but only against a background of consistent fulfillment, so the violation reads as
exceptional.

**B10 [GS] — the diagnostic use.** Applying B1–B7 to one's own prose *exposes missing
science*: the gaps in the argument surface as places where no old information links a
sentence to its predecessor. When revision cannot proceed without supplying a connection,
that connection is a claim the paper had failed to make. Treat every such point as a
substantive finding, not a wording problem.

### Diagnostic procedure (run on every Results paragraph)
1. List the topic position of each sentence. Do they name one continuing story?
2. List the stress position of each sentence. Is each one the thing worth emphasizing?
3. List the verbs. Do they name actions, or are they all copulas?
4. Any sentence where 1–3 fail is either badly built or is hiding a missing connection.

---

## C. Points of style (Whitesides)

**C1.** Do not use nouns as adjectives. "formation of ATP", not "ATP formation".
**C2.** "This" must always be followed by a noun, so its reference is explicit.
"This observation leads us to conclude", not "This leads us to conclude".
**C3.** Describe experimental results uniformly in the past tense.
**C4.** Use the active voice whenever possible — subject to B3 (topic position wins).
**C5.** Complete all comparisons. "higher using bromine **than chlorine**", not "higher
using bromine".

---

## D. Figures

**D1 [TSR-F].** Know your message; know your audience. One figure answers one question.
**D2 [TSR-F].** Do not trust the defaults; use color deliberately and accessibly; do not
mislead with axis choices.
**D3 [W]/[A2].** A figure is the primary carrier of a result. If a table and a figure in the
same section present the same trend, delete the table.
**D4.** Every figure caption states, in its first sentence, the claim the figure supports —
not a description of what is plotted.

---

## E. Project-specific editorial positions (locked earlier, retained)

**E1.** Report the threshold **plateau** (t = 1.0–2.5), not the peak. Canonical point t = 2.5.
**E2.** Quote the **rank correlation** for the tier gradient (ρ = −0.400, p = 0.0044), not
the binary odds ratio. OR = 6.48 is superseded and must not appear.
**E3.** The **cluster bootstrap** (resampling whole systems) is the primary uncertainty
estimate, not the per-variant bootstrap.
**E4.** State the AlphaMissense AUC gap **up front, as the premise** of the orthogonality
argument — never bury or hedge it.
**E5.** A negative structural call is uninformative about pathogenicity. Frame every
application as **enrichment of a candidate list, never exclusion** of a variant.
**E6.** Mechanism (how) and strength (how much) are reported separately from any clinical
pathogenicity verdict. Structural disruption and pathogenicity are orthogonal — a result
demonstrated, not assumed.
**E7.** Interpretation belongs in the Results of *this* paper. It is a methods/benchmark
paper; the interpretation of each benchmark result is the contribution, not a violation of
Results/Discussion separation.
**E8.** Report nulls at full strength. The tier ablation delta is exactly 0.0000 and is
stated as such.

---

## F. Standing numerical rules

**F0 — the denominator rule.** Every rate in the paper states its denominator, and the
denominator must match the question. Detection/firing rates are *mechanism* questions and
run over all variants of the expected mechanism class, **VUS included** — never filtered to
pathogenic/benign. Precision and pathogenicity-gradient statistics are *clinical-label*
questions and are necessarily restricted to labeled variants. Mixing the two produced a
spurious "fold arm 3/7" (correct value 6/10) and a spurious "sensitivity 0.52." Whenever a
rate appears, name the population it is computed over in the same sentence.

**F1.** Never transcribe a number from memory or from a previous draft. Every number in the
manuscript is regenerated from the canonical scored table by the verification script.
**F2.** The per-axis ANNOTATED counts (15/11/10) and GRADEABLE counts (13/10/10) are both
correct under different definitions. They are never reconciled to each other. Every use
states which definition is in force.
**F3.** Fig. 3c's monomer-fold arm (n = 14; 3 systems; 10 destabilizers / 4 intact) is a
direction check and is a different metric from the 21/27 structural-agreement figure in
Fig. 3b. Wherever both appear, one sentence states the difference explicitly.
**F4.** BRCT pathogenicity is read from `supplement/brct/brct_foldx_concordance.csv`
column `clinvar_germline`. The canonical `phenotype` field is NULL for
`system == "brca1_brct"` and silently produces wrong subcounts.
