# MAVIS v7 — Canonical Benchmark Ledger (post-verification, post-relaxed-tier, post-BRCT-expansion)

> **v6 update (BRCT monomer-fold expansion).** The benchmark has been extended from 44→56 variants
> (11→12 systems) by adding the **BRCA1 tandem-BRCT** monomer-fold system (12 variants, PDB 1JNX),
> grounded in **measured** GdmCl unfolding ΔΔG_U–F (Rowling, Cook & Itzhaki 2010, JBC 285:20080).
> This expands the previously-thin monomer-fold destabilizer arm from **n=2 → n=10** and adds a
> quantitative FoldX-vs-measured anchor. §§1–9 below describe the frozen 44-variant core (unchanged);
> **§10 records the v6 expansion and recompute.** Canonical variant table: `benchmark_variants_v6.csv`.

**Purpose.** The single source of truth for what the 44-variant benchmark *is* and why each call was
made. Supersedes the scattered append blocks for reference purposes (the append blocks remain the
detailed derivation record). For a methods paper where the benchmark is the contribution, this
ledger is a core deliverable, not hygiene.

## 1. Composition
- **44 variants / 11 PPI systems.** Roles (post mech-control retirement): benign 10, pathogenic 26,
  pathogenic_gof 6, pathogenic_lof 2.
- **No `mechanism_control` category** — the 11 originally-mislabeled rows are *mislabel corrected*
  (not reassigned): KRAS G12D/G12V/Q61H + PIK3CA H1047R → `pathogenic_gof`; BRCA1 C61G, MLH1 R755W,
  MSH2 G674R/C697F, VHL W117R/Y98H, HBB E6V → `pathogenic`.

## 2. Three ground-truth tiers (per axis: monomer, complex-fold, binding)
- **STRICT** (`expected_ddg_{axis}`): primary-literature-grounded. Token ∈ {stab, neutral, destab,
  unknown}. `unknown` = no groundable direct measurement (Definition B; FoldX self-prediction and
  bare structural-position inference are NOT admissible grounding).
- **RELAXED** (`expected_ddg_{axis}_relaxed`): strict + directional promotions of contracted axes
  that have indirect-but-specific support (coupled readout | sibling analogy | off-axis measurement
  that speaks to THIS axis). Structural-position inference excluded. One tier.
- **FUNCTIONAL** (not yet a column): off-axis quantitative data that doesn't speak to any structural
  axis (e.g. CaM Ca²⁺-affinity) — noted in evidence, reserved for a possible future functional axis.

## 3. Verification impact (the before/after)
- Graded structural axes: **96 → 80** (net −16; gross 18 contracted to unknown, 2 recovered to neutral).
- Disruptor tokens (destab/stab/mild_destab): **39 → 26** (a third removed as ungroundable or
  wrong-direction).
- Roles reassigned: **11** (mech-control retirement).
- **structural_agreement essentially unchanged (thresholded 0.77 @ t=2.5; the directional 0.773 is superseded — not reproducible from released columns, see §6):** the benchmark was corrected
  without moving the headline — the prior number was right for partly-wrong reasons, now right for
  right reasons.

## 4. The 18 contracted axes — relaxed disposition (6 promoted / 12 stay-unknown)
| variant | axis | strict | relaxed | basis |
|---------|------|--------|---------|-------|
| MSH2 A636P | mono | unknown | **destab** | COUPLED: instability/reduced-expression (Ollila) |
| MSH2 A636P | fold | unknown | unknown | no independent directional readout |
| MSH2 C697F | mono | unknown | unknown | no directional readout |
| MSH2 C697F | fold | unknown | unknown | no directional readout |
| MSH2 N127S | mono | unknown | **neutral** | COUPLED: MMR-proficient functional control |
| MSH2 N127S | fold | unknown | **neutral** | COUPLED: preserved assembly *(unevaluable: no partner)* |
| MSH2 G322D | mono | unknown | **neutral** | COUPLED: near-neutral polymorphism |
| MSH2 G322D | fold | unknown | **neutral** | COUPLED: preserved assembly *(unevaluable: no partner)* |
| VHL W117R | mono | unknown | unknown | structural-position only (excluded) |
| VHL W117R | fold | unknown | unknown | structural-position only (excluded) |
| VHL W117R | bind | unknown | unknown | not in Kishida panel; HIF-axis lesion |
| VHL Y98H | mono | unknown | **neutral** | COUPLED: retained VBC assembly proxy |
| TNNI3 R162W | bind | unknown | unknown | inference-only (Zhou 2013 not in hand) |
| CaM D96V | mono | unknown | unknown | off-axis Ca²⁺ ≠ fold (NMR: fold intact) |
| CaM N98S | mono | unknown | unknown | off-axis Ca²⁺ ≠ fold |
| CaM F142L | mono | unknown | unknown | off-axis Ca²⁺ ≠ fold |
| SMAD4 D351H | fold | unknown | unknown | Shi = oligomerization (binding axis), not fold |
| SMAD4 R361C | fold | unknown | unknown | Shi = oligomerization, not fold |

Promotions: 6 (1 destab, 5 neutral), all COUPLED-class. Conservative by design.

## 5. pLDDT interface-exclusions (4 variants — apply to BOTH metrics, BOTH pipelines)
| variant | partner | pLDDT | excluded axis | batch |
|---------|---------|-------|---------------|-------|
| TNNI3 R145G | tnnc1 | 54.93 | fold+binding | B9 |
| TNNI3 R145Q | tnnc1 | 54.93 | fold+binding | B9 |
| CaM D96V | cacna1c | 69.45 | fold+binding | B10 |
| SMAD4 D351H | smad3 | 57.37 | binding | B11 |

## 6. Locked headline metrics (Pipeline 1, t=2.5, pLDDT-reconciled)
- **structural_agreement (strict, thresholded sweep):** 0.77 @ t2.5 (range 0.76–0.80).
- **structural_agreement (directional, 1.0 floor) — SUPERSEDED:** the ledger’s strict 0.773 (51/66) / relaxed 0.757 (53/70) are **not reproducible from the released `mavis_v7_concordance_annotated.csv`** (the denominator counts axes per-partner for multi-chain systems plus the §5 pLDDT exclusions; reconstructions bracket 0.79–0.82). Report the thresholded 0.77 as the primary `structural_agreement`; the directional is omitted from the release unless regenerated from the grading script with its axis convention documented.
- **mech_consistency:** 0.70 raw → **0.73 pLDDT-reconciled** @ t2.5.
- **Pipeline 2:** graded locally — mech_consistency identical to P1 (both 0.698 raw), so the neighborhood pipeline adds no mechanism resolution; the neighborhood tier degrades the pathogenicity gradient (tier OR 6.48 → 4.00; neighborhood-elevated subset anti-enriched, OR 0.26) → tested-and-rejected alternative. (The earlier ledger figure “OR 0.48” does not reproduce under either natural 2×2 and is superseded by these.)

## 7. Per-system one-liners (mechanism character)
- **MLH1-PMS2** (Kosinski): mixed; L749P fold→unknown corrected.
- **HBB tetramer/dimer** (Kiger): W37 series R-state contacts — FoldX FN against T-state 2HHB.
- **KRAS-CRAF** (Hunter, closed): G12D/G12V GoF, low-RAF-affinity; Q61H direction was uncertain.
- **BRCA1-BARD1**: C61G (RING/Zn) mislabel-corrected to pathogenic.
- **PIK3CA-PIK3R1**: E545K/E542K released-autoinhibition GoF → FoldX reads as stabilizing.
- **VHL-ElonginC**: L158Q true interface disruptor; W117R fully contracted (HIF-axis); Y98H clean
  ElonginC-silent control.
- **MSH2-MSH6**: regulatory/DNA-binding system, zero true partner disruptors; A636P Walker-A.
- **TNNI3-TNNC1**: regulatory system, zero partner disruptors; R145G/Q binding preserved.
- **CaM-Cav1.2**: Ca²⁺-sensing system; CDI loss with interface preserved-to-strengthened; FoldX
  correctly "stabilizing/neutral" yet blind to the sensing lesion.
- **SMAD4-SMAD3**: R361C clean interface TP; D351H interface FN (AF2 mis-modeled, pLDDT-excluded);
  I500T/V structurally-silent post-translational GoF.

## 8. Cross-cutting principles confirmed by the campaign
- Three-axis separation is load-bearing; "right for the wrong reasons" is a first-class failure mode.
- `unknown` is a valid value; primary literature required; FoldX self-prediction and structural-
  position inference are not admissible grounding (strict tier).
- Off-axis evidence promotes a relaxed token ONLY if it speaks to that axis's direction.
- pLDDT gating applied uniformly across structural_agreement and mech_consistency, both pipelines.

## 9. Artifacts (this session)
- `mavis_v7_results_corrected_v2.csv` — strict tokens, post-recompute.
- `mavis_v7_results_relaxed.csv` — adds parallel `expected_ddg_*_relaxed` columns.
- `mavis_v7_concordance_annotated.csv` (released canonical name; formerly `mavis_v7_concordance_v5_reconciled.csv` — confirmed the same artifact: the reconciled file is absent from disk and the annotated file carries the full raw + pLDDT-reconciled mech_consistency columns across all thresholds) — P1 mech_consistency raw + pLDDT-reconciled, all thresholds.
- `batch12_recompute_bundle.py`, `mech_consistency_plddt_patch.py`, `relaxed_regrounding_walk.py` —
  the auditable edit/grade/promotion logic.
- Reports: Batch-12 cross-system, Track-B mech_consistency, reconciliation summary, P2 runbook,
  this ledger.

## 10. v6 expansion — BRCA1 tandem-BRCT monomer-fold system (44→56 variants, 11→12 systems)

**Motivation.** The monomer-fold destabilizer arm was the one under-supported axis (n=2: BRCA1 C61G,
MLH1 R755W). BRCT adds a graded ladder with **direct measured** subunit stability under Definition B.

**System.** `brca1_brct` — PDB 1JNX, tandem BRCT (residues 1646–1863), chain X, BRCA1 native numbering.
Ground truth = ΔΔG_U–F from Rowling, Cook & Itzhaki 2010 J Biol Chem 285:20080 (PMC2888420, Table S1;
GdmCl equilibrium unfolding, WT ΔG_unf = 10.56 kcal/mol at 10 °C). Truth tokens assigned from measured
ΔΔG via the pipeline's own `discretize_ddg` (t=1.0), held fixed; FoldX ΔΔG monomer is the prediction.

**12 variants.**
- Destabilizers (measured ΔΔG > 1.0, → strict monomer/fold `destab`): Y1853C (6.04), A1843P (4.89),
  V1736A (4.20), M1783T (3.73), V1808A (2.40), V1665M (2.22), R1751Q (1.57, benign), L1664P (1.18, benign).
- Neutral controls: M1663K (−0.03), P1806A (0.06).
- Fold-intact / function-lost pathogenics (pSer groove): R1699L (−0.99), R1699Q (−1.83).
Thesis illustration: strong destab → pathogenic; mild *real* destab → benign (R1751Q, L1664P);
fold-intact groove lesion → pathogenic (R1699L/Q) — i.e. disruption ⊥ pathogenicity, within one domain.

**Monomer-fold destabilizer arm: n=2 → n=10** across 3 systems (brca1_bard1 C61G, mlh1_pms2 R755W,
+ 8 BRCT destabilizers).

**FoldX vs measured (quantitative anchor).** MAVIS FoldX 5.1 (RepairPDB → BuildModel n=5) on 1JNX:
- Spearman ρ = **0.72** (fold axis, n=10, p=0.019); ρ = 0.67 over all 12.
- My FoldX vs Rowling's own reported FoldX: Pearson r = **0.83** (n=11) — confirms correct FoldX operation.
- BRCT fold-axis categorical accuracy: **10/12 = 0.833** (t=1.0/1.5/2.0).
- **V1736A** is a documented FoldX under-prediction (FoldX 1.24 vs measured 4.20; reproduces Rowling's
  own FoldX 0.02 — buried Val→Ala cavity). It fires correctly at the primary t=1.0 but is a
  false-negative at the stricter t=2.5; report both. This is the case that motivates measurement-anchoring.

**Recompute — structural_agreement (repo `evaluation.level3_mechanism_axis`, faithfully reproduced).**
Engine first reproduced the frozen 44-set headline exactly (85/109 = 0.780 @ t=1.0; 88/109 = 0.807 @ t=2.5).

| threshold | 44-set (frozen) | 56-set (with BRCT) | fold axis (56) |
|-----------|-----------------|--------------------|----------------|
| 1.0 (production) | 85/109 = 0.780 | **95/121 = 0.785** | 24/35 |
| 1.5 | 86/109 = 0.789 | 96/121 = 0.793 | 25/35 |
| 2.0 | 88/109 = 0.807 | 98/121 = 0.810 | 25/35 |
| 2.5 (ledger convention) | 88/109 = 0.807 | 98/121 = 0.810 | 25/35 |

**Headline is stable** (0.780 → 0.785 @ t=1.0): 12 measurement-grounded variants do not dilute the
metric — they nearly double the fold-axis evidence base (23 → 35 evaluable). The two BRCT fold-axis
misses at t=1.0 are both FoldX *over*-predictions: P1806A (neutral control, FoldX 1.67) and R1699Q
(fold-intact groove pathogenic, FoldX 1.77 vs measured −1.83) — the latter exactly the failure
measurement-anchoring is designed to catch.

> Absolute `structural_agreement` remains reproduction-fragile (§6): the delta 44→56 is computed from a
> single self-consistent engine run, so the +0.005 is the trustworthy quantity, not a change in the
> canonical 0.77 figure.

**v6 relaxed tier confirmed unchanged.** The BRCT expansion does not alter the §4 strict/relaxed
dispositions; `benchmark_variants_v6.csv` carries the §4 tokens verbatim (verified token-for-token) plus
the new BRCT rows.

## 11. Artifacts (v6 / BRCT session)
- `benchmark_variants_v6.csv` — 56 variants / 12 systems, strict + relaxed tiers + FoldX ddg_monomer column (canonical).
- `scored_56var_with_brct.csv` — full scored dataframe (checkpoint) driving the recompute.
- `brct_foldx_ddg.csv` — FoldX 5.1 output, 12 BRCT variants (5 runs each).
- `brct_foldx_concordance.{csv,md}` — FoldX-vs-measured concordance + interpretation.
- `brct_foldx_vs_measured.png` — measured vs predicted ΔΔG scatter (V1736A annotated).
- `v6_recompute_results.md` — full recompute record (this §10).
- FoldX run package: `1JNX_processed.pdb`, `individual_list.txt`, `run_brct_foldx.py`.
