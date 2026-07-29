# COMAVI v7 — Canonical Benchmark Ledger (post-verification, post-relaxed-tier, post-BRCT-expansion, post-interface-expansion)

> **v7 update (fold-neutral interface disruptors).** The benchmark has been extended from 56→61 variants
> (12→14 systems) by adding two interface systems whose pathogenic variants disrupt a protein–protein
> interface with the **monomer fold intact** — **VWF A1–GPIbα** (PDB 1SQ0) and **CFH FH1-4–C3b** (PDB 2WII) —
> plus matched retained-binding and benign controls. This supplies the first within-benchmark proof point
> for the interface-resolution claim. Canonical 61-set headline (t=2.5): mechanism-consistency **0.72**,
> structural agreement **0.77** (92/120). **§12–13 record the v7 expansion and recompute; §15 records the
> v7.1 grading-rubric correction that produced these values.**
> Canonical dataset: `scored_61var_canonical.csv`.
>
> **v6 update (BRCT monomer-fold expansion).** The benchmark was extended from 44→56 variants
> (11→12 systems) by adding the **BRCA1 tandem-BRCT** monomer-fold system (12 variants, PDB 1JNX),
> grounded in **measured** GdmCl unfolding ΔΔG_U–F (Rowling, Cook & Itzhaki 2010, JBC 285:20080).
> This expands the previously-thin monomer-fold destabilizer arm from **n=2 → n=10** and adds a
> quantitative FoldX-vs-measured anchor. §§1–9 below describe the frozen 44-variant core (unchanged);
> **§10 records the v6 expansion and recompute.** v6 variant table: `benchmark_variants_v6.csv`.

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
- **structural_agreement (directional, 1.0 floor) — SUPERSEDED:** the ledger’s strict 0.773 (51/66) / relaxed 0.757 (53/70) are **not reproducible from the released `comavi_v7_concordance_annotated.csv`** (the denominator counts axes per-partner for multi-chain systems plus the §5 pLDDT exclusions; reconstructions bracket 0.79–0.82). Report the thresholded 0.77 as the primary `structural_agreement`; the directional is omitted from the release unless regenerated from the grading script with its axis convention documented.
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
- `comavi_v7_results_corrected_v2.csv` — strict tokens, post-recompute.
- `comavi_v7_results_relaxed.csv` — adds parallel `expected_ddg_*_relaxed` columns.
- `comavi_v7_concordance_annotated.csv` (released canonical name; formerly `comavi_v7_concordance_v5_reconciled.csv` — confirmed the same artifact: the reconciled file is absent from disk and the annotated file carries the full raw + pLDDT-reconciled mech_consistency columns across all thresholds) — P1 mech_consistency raw + pLDDT-reconciled, all thresholds.
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

**FoldX vs measured (quantitative anchor).** COMAVI FoldX 5.1 (RepairPDB → BuildModel n=5) on 1JNX:
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

## 12. v7 expansion — fold-neutral interface disruptors (56→61 variants, 12→14 systems)

**Motivation.** The v6 benchmark could demonstrate the binding axis only on variants that *co-destabilize
the fold* (all 6 v6 binding-axis hits had median Grantham 113 with concurrent fold ΔΔG). This left the
central interface claim — that COMAVI resolves binding-interface disruption where a monomer-only tool
cannot — as an architecture argument without a within-benchmark proof point. v7 adds two systems whose
pathogenic variants disrupt a protein–protein interface **with the monomer fold left intact**, plus
matched retained-binding / benign controls.

**Two systems (+5 variants, all with a binding-axis structure).**
- **VWF A1–GPIbα** — PDB **1SQ0** (von Willebrand factor A1 domain × platelet GPIbα).
  Numbering offset **VWF prepro −763** (mature = prepro − 763). K_D ground truth from
  Tischer / Moon-Tasson / Auton 2025 *J Thromb Haemost* 23(4):1215–1228, **PMID 39756657**.
  - **R1334Q** (pathogenic type-2M): mono **0.0051**, fold **4.7285**, bind **4.1078** — the
    signature pure-interface disruptor (fold-neutral monomer, fired binding axis).
  - **A1381T** (benign polymorphism): mono **−0.0984**, fold **−0.1061**, bind **−0.0015** — flat on
    all three axes (specificity control).
- **CFH FH1-4–C3b** — PDB **2WII** (complement factor H domains 1–4 × C3b).
  Numbering offset **CFH prepro −18**. K_D ground truth from Pechtl 2011 *J Biol Chem* 286(13):11082–90,
  **PMID 21270465**.
  - **R78G** (pathogenic aHUS): mono **0.9297**, fold **4.8753**, bind **6.1211**; measured K_D **>35 µM** —
    pure-interface disruptor, fold-neutral by the monomer axis, binding axis fired.
  - **R53H** (pathogenic aHUS, retained binding): mono **0.4688**, fold **1.3156**, bind **0.0075**;
    measured K_D **~12 µM** — pathogenic but binding-competent → correct binding-axis *negative*.
  - **I62V** (benign polymorphism): mono **0.7140**, fold **0.3759**, bind **0.0052**; measured
    K_D **10–14 µM** — benign, binding-competent control.

**Structural deviation (documented for methods).** 2WII is 14,081 atoms and stalled RepairPDB under
host memory pressure; C3b was trimmed to residues within **15 Å of chain C** (FH kept whole),
5,144 atoms (37%), running in ~9 min. AnalyseComplex on the partner split across ≥2 chains sums all
interface pairs. FoldX 5.1, COMAVI-faithful recipe (RepairPDB → BuildModel n=5 fold ΔΔG →
AnalyseComplex WT+5 mut binding ΔΔG; DDG_DESTAB=1.0).

**Named exclusion — applicability domain (5 VWF type-2M variants).** Five additional VWF type-2M
candidates — **L1282R, D1302G, V1360A, V1439M, I1425T** — were excluded *before scoring*, from
mechanism + ground truth (not from any COMAVI output): their loss-of-function is conformational /
mechanotransductive (shear-gated A1 exposure), a mode not reachable from a static A1–GPIbα interface
model. This is the benchmark's one principled exclusion; structurally-silent out-of-scope variants are
**retained** as specificity controls.

**61-set recompute (canonical basis: frozen 44-annotated sweep + 5-new operative-ΔΔG sweep, reconciled
with the tier axis).** Headline at the canonical calling threshold t=2.5:
- **mechanism-consistency 0.72** (0.7234; graded n=47), **structural agreement 0.77** (92/120 = 0.7667).
  *(Superseded by §15: the pre-correction values were 0.71 / 0.7083 at graded n=48 and 0.76 / 92/121.)*
- Threshold sweep (MC / SA): t1.0 0.553 / 0.725 (87/120); t1.5 0.660 / 0.775 (93/120);
  t2.0 0.681 / 0.767 (92/120); **t2.5 0.723 / 0.767 (92/120)**; tSAP 0.734 / 0.758 (91/120).

**Four-way structural-agreement decomposition (t=2.5, canonical basis):**
monomer **14/16** + complex-fold **20/25** + binding **24/32** + tier **34/47** = **92/120 = 0.767**.
*(Pre-§15 the tier term was 34/48, giving 92/121 = 0.760; only the tier denominator moved.)*

**Tier gradient (authoritative, 61-set):** T1 **14/14 (100%)**, T2 **13/18 (72%)**, T3 **7/10 (70%)**,
T4 **3/7 (43%)**; Fisher **OR = 3.78, p = 0.080**; Spearman **ρ = −0.400, p = 0.0044**
(strong tiers 1–2: 27/32; weak tiers 3–4: 10/17). The rank-correlation gradient survives the
expansion while the binary strong-vs-weak dichotomy loses Fisher significance — the ρ framing is
reported as primary, consistent with the tier being an orthogonal structural-disruption evidence
axis rather than a pathogenicity classifier.

**Physical validation unchanged** (independent of the expansion): BRCT fold ρ = 0.72 (n=10, p=0.019);
Hb binding ρ = 0.90 (n=7, p=0.037).

## 13. Artifacts (v7 / interface-disruptor session)
- `scored_61var_canonical.csv` — 61 variants / 14 systems, canonical basis (released canonical dataset).
- `vwf_axes_results.csv`, `cfh_axes_results.csv` — FoldX 5.1 three-axis output, new-5 variants.
- FoldX run packages: `1SQ0` (VWF A1–GPIbα), `2WII` trimmed (CFH FH1-4–C3b, 15 Å of chain C).
- Regenerated main figures on the 61-set: F2 (headline SA/MC sweep + tier gradient),
  F3 (per-axis competency + mechanism-class), F5 (AlphaMissense vs tier).

## 14. Per-axis count definitions — reconciliation (authoritative)

Three different per-axis triples circulate in the project documents and drafts. **All three are
arithmetically correct; they count different things.** This section fixes the definitions and gives
the canonical values, each reproduced by loading `scripts/apply_concordance_v5.py` and calling its
own `_operative_axis` / `_axis_check` on `reference_outputs/scored_61var_canonical.csv` (t = 2.5).

| # | Quantity | mono | cx-fold | binding | Definition |
|---|---|---|---|---|---|
| A | **Annotated destabilizers** | **10** | **19** | **15** | Rows whose `expected_ddg_*` ground-truth token is `destab`. No prediction-side gating. This is a property of the *curated benchmark*, not of any pipeline run. |
| B | **Gradeable axes** | **16** | **25** | **32** | Axes that enter the structural-agreement denominator: prediction present, axis-evaluability gate passed, and the prediction's internal CI95 excludes 0. Includes neutral and stabilizing controls. Sums with the tier axis (48) to the headline denominator **121**. |
| C | **Graded disruptors** | **2** | **10** | **13** | Intersection of A and B — annotated destabilizers that are *also* gradeable at t = 2.5. |

### Which to quote where

- **§2.2 / Methods, describing benchmark composition** → **A (10 / 19 / 15)**. The claim is "the
  curated set contains this many disruptors per axis," which is threshold-independent.
- **§3.3 / per-axis agreement table** → **B (16 / 25 / 32)**, already correct in the manuscript, with
  the existing prose explaining that B exceeds C because the denominator includes controls.
- **C is not a headline quantity** and should not appear as "the axis is supported by n disruptors" —
  it is threshold-dependent and, on the monomer axis, misleadingly small for a reason that is an
  artifact of file layout, not of the science (below).

### The monomer-axis trap

Column `ddg_monomer_distinguishable_internal_from_0` is **NaN for all 12 BRCT rows** in the canonical
CSV: BRCT monomer ΔΔG was computed in the v6 supplement run
(`supplement/brct/brct_foldx_ddg.csv`, n = 5 replicates with per-run values and SD) and merged in
without that derived CI column. Two consequences:

1. **`bool(NaN)` is `True` in Python.** A naive gate check written as `bool(row.get(col, False))`
   silently *admits* all 12 BRCT rows, giving monomer gradeable = 28 and disruptors = 10. The shipped
   `compute_structural_agreement` reads the value with a `False` default, so a genuinely missing
   column is excluded and monomer gradeable = **16** — which is what the manuscript table reports and
   what the headline 92/120 decomposition uses (§15). Any reimplementation of the gating must reproduce the
   `False`-default semantics; a plain truthiness test does not.
2. **Definition C therefore reads 2 on the monomer axis**, i.e. the pre-BRCT count — not because the
   BRCT expansion failed, but because the CI column those rows would need is absent from the canonical
   file. The expansion's evidentiary weight is carried by the measured-ΔΔG anchor (ρ = 0.72, n = 10
   core/fold sites, p = 0.019) and the 0.83 categorical accuracy, both computed in the supplement,
   not by C. **This is the reason C must not be quoted as arm support.**

Optional cleanup (not required for the paper, and it would touch a frozen number): backfill
`ddg_monomer_distinguishable_internal_from_0` for the BRCT rows from the per-replicate SDs already in
`brct_foldx_ddg.csv`, which would move monomer gradeable 16 → 28 and shift the headline denominator.
**Do not do this before submission** — the frozen 92/120 = 0.767 (§15) is computed on the current file.

### BRCT destabilizer-count double-threshold (verified, not a defect)

`supplement/brct/brct_foldx_concordance.csv` has `measured_destab` **True for 4** rows while the
canonical CSV annotates **8** BRCT rows `destab`. Both are correct at their own documented cut:

- `measured_destab` uses the **2.5 kcal/mol** cut of the standalone classification table
  (sensitivity/specificity 0.75, TP 3 / FN 1 / FP 2 / TN 6) — 4 of 12 measured ΔΔG_U–F values exceed 2.5.
- `expected_ddg_monomer` uses the pipeline's **primary 1.0 kcal/mol** operating point — 8 of 12 exceed 1.0.

Verified: applying `destab if measured > 1.0, stab if < −1.0, else neutral` reproduces all 12
canonical `expected_ddg_monomer` tokens **exactly** (8 destab / 3 neutral / 1 stab, no mismatches).
The four rows that differ (V1808A 2.40, V1665M 2.22, R1751Q 1.57, L1664P 1.18) are the
threshold-borderline band already discussed in `brct_foldx_concordance.md`.

## 15. v7.1 grading-rubric correction — `structurally_uncommitted` fall-through (authoritative)

**The bug.** `derive_expected_mech_class` had no branch for a variant whose per-axis strict ground
truth is *entirely* `unknown` **and** whose topology is `away_from_interface`. Such a variant fell
through to the terminal `return "structurally_silent"`. That is a substantive mislabel, not a naming
quibble: `structurally_silent` is an affirmative claim that the axes were measured and found
inactive, and it is graded as such — the variant is expected to produce no structural signal, and
any signal COMAVI reports is scored as a false positive. A variant with no admissible evidence on
any axis supports no such expectation and must not be graded at all.

**The fix.** A new terminal class `structurally_uncommitted` is returned when no axis carries a
committed token and topology is not interface-positive. The class is excluded from
mechanism-consistency grading (like the pre-existing `interface_uncommitted_magnitude`) and its
tier axis is excluded from `structural_agreement`, because the tier axis in that metric is
conditioned on `expected_mech_class`.

**Scope on the canonical 61-set: exactly one variant, VHL W117R.** Its strict ground truth is
`unknown` on all three ΔΔG axes (`evidence_axes = none`); the ledger §4 disposition records its
monomer axis as *"structural-position only (excluded)"* — position inference is explicitly not
admissible grounding under §2. It was previously graded `inconsistent` at all five thresholds:
stable, and uniformly wrong for a reason that was an artifact of the fall-through.

**Effect on the headlines (t=2.5).**

| metric | before | after |
|---|---|---|
| mechanism-consistency | 0.7083 (graded n=48) | **0.7234** (graded n=47) |
| structural agreement | 92/121 = 0.7603 | **92/120 = 0.7667** |

Full sweep after correction (MC / SA): t1.0 0.553 / 0.725 · t1.5 0.660 / 0.775 · t2.0 0.681 / 0.767
· t2.5 **0.723 / 0.767** · tSAP 0.734 / 0.758. Bootstrap/Jeffreys 95% CIs on the corrected
canonical values: MC [0.606, 0.830]; SA [0.685, 0.835].

**What did NOT move.** The numerator is unchanged at every threshold — no prediction changed. The
per-axis decomposition confirms the change is confined to one tier axis: monomer 14/16, complex-fold
20/25, binding 24/32 are identical before and after; only the tier term moves 34/48 → 34/47. The
**tier pathogenicity gradient is untouched** (14/14, 13/18, 7/10, 3/7 = 100 / 72 / 70 / 43 %;
Spearman ρ = −0.400, p = 0.0044, n = 49) because the tier gradient carries no `expected_mech_class`
term. Figures 2, 3 and 5 regenerate from the canonical CSV without edits and self-report the
corrected values.

**Reproduction / validation of the correction.** Before applying the fix, the grading stage was
re-executed on the *unmodified* canonical table and reproduced all five stored `structural_agreement`
numerators and denominators **exactly** (87/121, 93/121, 92/121, 92/121, 91/121). Two conventions are
required to do so and are recorded here because they are easy to get wrong:
1. `discover_partners()` over the released wide table also returns the per-partner `_ci95_*` and
   `_distinguishable_*` columns as if they were partner chains; they must be filtered out.
2. Grading is applied to the **49 interaction rows only** (`expected_mech_class` non-null). The 12
   BRCA1-BRCT fold-expansion rows are deliberately excluded from interaction-axis grading; including
   them inflates the denominator to 133.

**Known gap — the verification harness does not cover this.** `verification/verify_stage6.py` passes
11/11 unchanged after the correction, and that pass is **not** evidence for the released numbers. The
harness runs against `inputs/intermediate/comavi_v7_results_with_nbhd.csv` — the 44-variant
generation — and its `EXPECTED_V5` constants date from that generation (SA t2.5 0.718, MC t2.5 0.761).
VHL W117R is present in that file, but with its **pre-re-curation** ground truth
(`expected_ddg_monomer/fold_complex = destab`, `binding = neutral`, `evidence_axes = both`), which
derives `fold_mechanism` and never reaches the fall-through. Diffing the 44 shared variants between
that intermediate and the canonical table shows **14 of 44 had their per-axis ground truth changed**
between generations — overwhelmingly downgrades to `unknown` under the tightened §2 standard, and
several (MSH2 A636P, TNNI3 R145G/R145Q/R162W, SMAD4 D351H/R361C, VHL W117R) change expected class.
That re-curation is intentional and documented (§2–§4), but its consequence is that the harness
validates the older, looser curation and is silent on the released canonical numbers.

**Resolved — Track C.** `verify_stage6.py` now takes `--canonical` and adds a third track that
verifies the released table directly (`EXPECTED_CANONICAL`). There is no 61-variant *pipeline input*
in the repo — the canonical table is itself the released product — so Track C cannot re-run the
pipeline end to end. Instead it (1) re-derives the grading columns from the released table using the
pipeline's own `derive_expected_mech_class` / `classify_axis_status` /
`grade_mechanism_consistency` / `compute_structural_agreement`, and asserts they reproduce the STORED
columns exactly at all five thresholds — which is what guards against the table drifting from the
code that produced it; (2) checks the pooled MC/SA headlines and the graded n; (3) asserts the tier
gradient is unchanged, since it must not depend on grading. Track B is retained as a regression test
on the 44-variant generation and is explicitly labelled as such. Full run: **33/33 pass**
(11 Track B + 22 Track C).

    python3 verification/verify_stage6.py       --intermediate inputs/intermediate/comavi_v7_results_with_nbhd.csv       --am inputs/AM_variants_comavi_mechanism_test.xlsx       --scripts-dir scripts       --canonical reference_outputs/scored_61var_canonical.csv

---

## 16. Benchmark stress tests (v7.1, authoritative)

Four robustness tests on the released canonical table, in
`verification/stress_tests.py`. All four re-derive expected mechanism class and
per-axis status from the ground-truth columns using the pipeline's own functions
in `apply_concordance_v5`, then grade against the stored predictions. FoldX is
never re-run; no ΔΔG is recomputed.

The scorer is gated: it must reproduce the stored per-variant grades exactly
(`mech_consistency_t25`) before any perturbation is applied. That assertion
fires on every run.

| Test | Result |
|---|---|
| **A. Permutation null** (n=2000, ground-truth block permuted as a unit) | MC observed **0.7234** vs null 0.4950 ± 0.0510, **p = 5×10⁻⁴** (floor for 2000 draws); SA observed **0.7667** vs null 0.5381 ± 0.0495, **p = 5×10⁻⁴** |
| **B. Leave-one-system-out** (14 gene systems) | MC range **0.7045–0.7619**; SA range **0.7522–0.7788**. No single system carries the result. |
| **C. FoldX replicate noise** (n=500, Gaussian at SE of the 5-replicate mean, mechanism labels **re-called** from perturbed energies) | MC **0.7256 ± 0.0162**; SA **0.7645 ± 0.0040**; mean 0.97 label flips per draw (max 4 of 49); **29 % of draws label-identical**; **41/49 variants never flip** |
| **D. Variant bootstrap** (n=2000) | MC 0.7234, 95 % CI **[0.606, 0.827]**; SA 0.7667, 95 % CI **[0.689, 0.839]** |

**Noise-fragile calls** (`verification_output/comavi_noise_fragility.csv`) — the
eight variants whose mechanism label is not fully stable under the force field's
own replicate scatter. Report these as a diagnostic, not as failures: they are
variants sitting near a calling threshold.

| Variant | Mechanism at t=2.5 | Flip rate |
|---|---|---|
| TNNI3 R145Q | Multimer fold destabilization at interface | 0.448 |
| PIK3CA H1047R | Structural variant — contact-driven (ΔΔG neutral) | 0.316 |
| CALM1 D96V | Interface variant (ΔΔG neutral) | 0.110 |
| PIK3R1 N564D | Multimer fold + PPI destabilization | 0.058 |
| MSH2 G674R | Both fold destabilization | 0.046 |
| MLH1 L749P | Multimer fold destabilization at interface | 0.034 |
| SMAD4 R361C | Both fold + PPI destabilization | 0.024 |
| TNNI3 R145G | Multimer fold destabilization at interface | 0.006 |

### Three design decisions worth preserving

1. **The permutation null permutes the whole ground-truth block as a unit**, not
   each column independently. Independent column shuffling fabricates mechanism
   profiles that never occur in nature (e.g. binding-lost with no interface
   contact) and yields a falsely low null.
2. **Test C must re-call the mechanism label from the perturbed energies**
   (`recall_mech=True`). Grading perturbed ΔΔG against the *stored* label makes
   mechanism-consistency mathematically unable to move — the first run reported
   sd = 0.0000, which was an artifact of that mistake, not stability.
3. **Panel (c) of the supplementary figure is a stem plot, not a histogram.**
   Mechanism-consistency on n=47 graded variants is discrete in steps of
   0.5/47 = 0.0106; only 8 values are attainable, and binning invents gaps.

### Reproduce

    python3 verification/stress_tests.py \
        --canonical reference_outputs/scored_61var_canonical.csv \
        --scripts-dir scripts --out-dir verification_output
    python3 figures/src/figureS_stress_tests.py

The figure script carries its own geometric layout assertion (no label overlap,
no spine collision, nothing clipped) and exits non-zero rather than writing a
figure that fails it. Seeds are fixed in the script, so both are deterministic.
