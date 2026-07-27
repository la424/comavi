# MAVIS — Pre-Publication Checkpoint

> ## ⚑ HISTORICAL DOCUMENT — superseded by the v7 / 61-variant canonical set
>
> **This checkpoint is dated 2026-07-01 and records the v6-era (44/56-variant) state of the
> project. It is retained for the audit trail. Do not cite any number, file path, or open task
> from this document as current.**
>
> Current authoritative sources:
> - **Benchmark numbers** → `docs/MAVIS_v7_canonical_benchmark_ledger.md` §12
> - **Canonical scored set** → `reference_outputs/scored_61var_canonical.csv` (61 variants, 14 systems)
> - **Figures** → `figures/README.md`
> - **Results interpretation** → `docs/MAVIS_results_synthesis.md` (top banner governs)
>
> **What changed since this checkpoint (read this before §2–§8):**
> 1. **The monomer-fold expansion named as the gating task in §7 is DONE** (v6, BRCA1 tandem-BRCT,
>    PDB 1JNX, 12 variants, Rowling 2010 measured GdmCl ΔΔG_U–F). §2's "n=2 graded destabilizers"
>    and §7's expansion task are both obsolete. See ledger §10.
> 2. **The canonical benchmark file changed.** §3 and §6.4 name `mavis_v7_concordance_annotated.csv`
>    (44 variants); the canonical set is now `reference_outputs/scored_61var_canonical.csv` (61).
> 3. **The headline numbers moved** (§4 banner below has the current values).
> 4. **§8's "in sync with GitHub" is point-in-time** and does not describe the working tree now.
>
> **Only remaining open item project-wide:** CHD per-axis evaluability gating (§7, production-only —
> no effect on the methods paper).

**Date:** 2026-07-01  **Scope:** both papers (methods/benchmark + CHD application), evaluated in parallel.

**Verdict (as of the checkpoint date):** the pipeline is in strong shape and reproducible at the
concordance level. One CHD takeaway moved with the AlphaMissense parse fix (§5); benchmark takeaways
are unchanged. Two items remained before submission at that time — one substantive (the monomer-fold
cohort, **since resolved in v6**), one production-only.

## 1. Inputs — complete
_(Paths updated to the current repo layout; the v1.2 reorg moved the benchmark input to `inputs/raw/` and the CHD inputs to `inputs/chd/`. Content otherwise as of the checkpoint date.)_
- **Benchmark:** `inputs/raw/benchmark_variants_v5.csv` (44 variants / 11 PPI systems).
- **CHD:** `inputs/chd/chd_input_final.csv`, `inputs/chd/chd_input_per_system.csv` (bidirectional per-system),
  `data/MAVIS_CHD_variant_domain_interface_map.csv`.
- **Annotations:** `inputs/chd/variants_with_alphamissense_and_franklin_expanded.csv` (AlphaMissense + Franklin,
  transposition-corrected), `inputs/AM_variants_mavis_mechanism_test.xlsx`.
- **Structures:** AlphaFold-3 Server (multimers); experimental PDBs where AF2 confidence is insufficient
  (2HHB / HBB, 6XI7 / KRAS-RAF1, 1JM7 / BRCA1-BARD1). `structures/` is gitignored (obtain per docs).

## 2. Pipeline — works as intended, with known limits
The `mavis_v7` engine is byte-identical across the working and release trees; the concordance engine
round-trips byte-for-byte on the canonical CHD file; `archive/derivation/apply_chd_concordance.py` is idempotent on
already-correct input.

Remaining limits:
- ~~**Monomer-fold axis undersupported — n=2 graded destabilizers.** The methods paper's one substantive
  gating task (see §7).~~ **[RESOLVED in v6]** The BRCA1 tandem-BRCT expansion (PDB 1JNX, 12 variants,
  measured GdmCl ΔΔG_U–F from Rowling 2010) closed this arm. On the canonical 61-set the monomer-fold
  axis carries **10 annotated destabilizers across 3 systems** (brca1_brct 8, brca1_bard1 1,
  mlh1_pms2 1) and is anchored to measured free energies (Spearman ρ = 0.72, n = 10 core/fold sites,
  p = 0.019). See ledger §10.
- **Full FoldX-level reproducibility gap.** The pLDDT-reconciled structural layer was not persisted, so
  outputs regenerate from the locked structural CSV (concordance level), not from raw FoldX. State in methods.
- **Directional `structural_agreement` is not reproducible from released columns** → report thresholded 0.77.

## 3. Outputs — inventory
| Artifact | Path | Notes |
|---|---|---|
| Benchmark comprehensive **[SUPERSEDED]** | `reference_outputs/mavis_v7_concordance_annotated.csv` | 44 variants, per-axis, 1255 cols — v6-era. **Canonical set is now `reference_outputs/scored_61var_canonical.csv` (61 variants, 14 systems).** |
| CHD comprehensive | `reference_outputs/chd_concordance_results_FIXED.csv` | 384 rows × 213 cols |
| CHD collapsed | `reference_outputs/chd_concordance_collapsed.csv` + `MAVIS_CHD_concordance_collapsed.xlsx` | 144 variants, one row each |
| Per-variant overview | `reference_outputs/MAVIS_results_summary.xlsx` | 11 sheets (benchmark + CHD per-variant, distributions, candidates, control recall, orthogonal cases) |
| Narrative summary | `reference_outputs/MAVIS_results_summary.docx` | prose writeup + matching tables |
| Benchmark takeaways | `docs/MAVIS_v7_canonical_benchmark_ledger.md` | locked metrics, per-system notes, principles |
| Two-paper synthesis | `docs/MAVIS_results_synthesis.md` | benchmark + CHD results, interpretation, verification ledger |
| Methods / decisions | `docs/methods_metrics_sketch.md`, `docs/design_decisions.md` | |

## 4. Locked headline results

> **⚑ SUPERSEDED — v7 / 61-variant canonical set.** The benchmark headline below is the v6-era
> 44/56-variant record. Current authoritative numbers (recomputed from
> `reference_outputs/scored_61var_canonical.csv`; full record in ledger §12): **61 variants** (49 PPI
> across 13 complexes + 12 BRCA1-BRCT; 14 systems); structural_agreement **0.76** (92/121 at t=2.5),
> mech_consistency **0.71** (graded n=48); tier gradient **100 / 72 / 70 / 43** (T1–T4), Fisher
> OR = 3.78 (p = 0.080), Spearman ρ = −0.40 (p = 0.0044) — the rank-correlation framing is preferred.
> Physical validation unchanged: BRCT fold ρ = 0.72 (n=10), Hb binding ρ = 0.90 (n=7). AlphaMissense:
> 47 PPI with a score (36 pathogenic / 11 benign). No tier×FoldX AUC-improvement claim is made
> (Spearman ρ = 0.60 partial-independence only; the "0.81→0.87" observation did not reproduce on 61).

**Benchmark** (P1, t=2.5, pLDDT-reconciled): structural_agreement **0.77** (thresholded — primary);
mech_consistency **0.73** reconciled / 0.70 raw; tier gradient **100 / 81 / 70 / 33** (T1–T4);
AlphaMissense accuracy on confident calls 37/41 = 90%. Claim arms: PPI-disruption **13**, complex-fold
**11**, monomer-fold **2**, silence-not-benign **19**. All 6 GoF variants reframed as structurally silent.
Pipeline 2 (neighborhood) equals P1 on mechanism and degrades the pathogenicity gradient → tested-and-rejected.

**CHD:** 144 variants (133 patient + 11 ZIC3 controls); 45 evaluable, 99 unevaluable (disordered in every
context). Controls: 8/9 structure, 7/9 ddG. **Patient candidates: 8** — KPNA6 I498T at 4/4 (all four
channels), seven at 3/4. ROCK2 T367M is the single clean MAVIS-vs-AM disagreement; DVL2 D441Y the cleanest
novel candidate; H286R correctly scored silent (trafficking / NLS mechanism).

## 5. Did the recent AlphaMissense fix change any takeaways?
- **Benchmark: no.** The fix was CHD-only (kpna1 / kpna6 / tcf7l1); no benchmark gene, metric, or claim moved.
- **CHD: one.** Correcting the whole-gene AlphaMissense score→class transposition restored KPNA6 I498T's
  AM hit, lifting it 3/4 → **4/4** (`concordance` metric) and adding KPNA6 K424N as an 8th candidate
  (2/4 → 3/4). I498T is now AM-concordant (previously logged as AM-discordant). The central CHD thesis is
  intact (ROCK2 T367M, DVL2 D441Y, structural disruption ≠ pathogenicity); under strict pathogenic-only
  Franklin the novel candidates still sit at 3/4. Reflected across the CSVs, the summary workbook, and both
  narrative docs.

## 6. Resolved at this checkpoint (consistency pass)
1. Report thresholded `structural_agreement` **0.77** as primary; directional 0.773 superseded (not
   reproducible from released columns; ledger §6 updated).
2. Silence-not-benign = **19** (verified): non-benign variants with no `destab` token on any structural
   axis; the "29" was the all-silent count, which also counted the 10 benign variants.
3. Pipeline 2 = tested-and-rejected; report tier OR 6.48 → 4.00 and elevated-subset OR 0.26; the ledger's
   "OR 0.48" is superseded (does not reproduce under either natural 2×2).
4. Canonical benchmark file = `mavis_v7_concordance_annotated.csv` (formerly `v5_reconciled`, absent from
   disk; confirmed the same artifact — carries the full raw + pLDDT-reconciled mech_consistency columns).

## 7. Open before submission
- ~~**[Substantive] Expand the monomer-fold destabilizer cohort** (2 → 6–8 graded). Leads BRCA1 V11G,
  MLH1 Q542L, SMAD4 R361C, TNNI3 R162W were confirmed non-gradeable; productive path = unstable-hemoglobin
  variants (HBB already on 2HHB) and/or a new monomer-fold system. Reopens the recompute; gating for the
  methods paper.~~ **[RESOLVED in v6]** — closed via BRCA1 tandem-BRCT (PDB 1JNX, 12 variants) rather than
  hemoglobin. The recompute was reopened and completed: the canonical set went 44 → 56 → **61 variants**.
  No substantive item gates the methods paper.
- **[Production — STILL OPEN] CHD per-axis evaluability gating refinement** — gate each axis on its own
  context's pLDDT (monomer-fold ← monomer; complex/binding ← complex-position). No effect on the current
  45/99 split; matters on the full production run. **This is the only open item project-wide, and it does
  not affect the methods paper.**

## 8. Repository state
_Point-in-time as of 2026-07-01; see `git status` for the live state._

`la424/mavis`, branch `main`, in sync with GitHub as of this checkpoint date. The FoldX binary and
AlphaFold structures are gitignored — users supply their own (see `README.md` and `docs/`).
