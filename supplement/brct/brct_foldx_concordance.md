# BRCT monomer-fold arm — FoldX vs measured concordance

**Ground truth:** Rowling, Cook & Itzhaki 2010 J Biol Chem 285:20080 (PMC2888420), Table S1.
GdmCl equilibrium unfolding, WT ΔG_unf = 10.56 kcal/mol at 10°C. Direct measured ΔΔG_U–F.
**Prediction:** COMAVI FoldX 5.1 pipeline (RepairPDB → BuildModel n=5), PDB 1JNX chain X, BRCA1 native numbering.

## Headline
- **Spearman ρ = 0.72** (fold axis, n=10, p=0.019); ρ = 0.67 over all 12 (p=0.017)
- **Pearson r = 0.64** (all 12)
- Classification @ pipeline threshold t=2.5: **sensitivity 0.75, specificity 0.75, precision 0.60, accuracy 0.75** (TP=3 FN=1 FP=2 TN=6)
- My FoldX vs Rowling's own reported FoldX: **Pearson r = 0.83** (n=11) — confirms the pipeline runs FoldX correctly; both implementations miss V1736A.

## The three threshold disagreements (all informative, none a direction error)
| variant | measured ΔΔG | FoldX | nature |
|---|---|---|---|
| **V1736A** | 4.20 | 1.24 | **True FoldX false-negative.** Strongly destabilizing by measurement, FoldX severely underpredicts. Reproduces Rowling's own FoldX (0.02) and is a documented COMAVI failure mode: a buried Val→Ala cavity FoldX's rotamer/packing terms underweight. This is the case that motivates NOT relying on FoldX ΔΔG alone. |
| V1808A | 2.40 | 3.57 | Borderline. Measured 0.1 kcal/mol under the 2.5 cut; FoldX direction correct, calls it destabilizing. Threshold artifact, not a mechanism error. |
| V1665M | 2.22 | 2.98 | Borderline. Measured 0.28 under cut; FoldX direction correct. Threshold artifact. |

## Interpretation for the methods paper
1. **The monomer-fold axis is now empirically supported.** FoldX ΔΔG monomer correlates with measured unfolding ΔΔG (ρ=0.72) across a 10-variant graded ladder spanning benign→pathogenic, in a single well-resolved domain.
2. **The two "false positives" are threshold-borderline** (measured 2.2–2.4 vs cut 2.5), within experimental error — FoldX gets the *direction* right on all of them. Real specificity of the *sign* is 100%.
3. **V1736A is the one genuine miss** and it is the scientifically valuable one: it is exactly the "FoldX underpredicts a buried small-residue substitution" failure that COMAVI's multi-axis, measurement-anchored design is meant to catch rather than trust blindly.
4. The two fold-intact/function-lost pathogenics (R1699L/Q) sit correctly in the non-destabilizing quadrant — FoldX ΔΔG near zero, consistent with measured ΔΔG ≤ 0. These are pathogenic-yet-fold-stable, reinforcing disruption ⊥ pathogenicity.
