# COMAVI main figures (methods/benchmark paper, 61-variant set)

<!-- COMAVI_ISDS_V1_FIGURE_BANNER -->

## Current publication figure suite

The current COMAVI figures are under [`isds_v1/`](isds_v1/), with alt text,
checksums, plotting source, and a regeneration command. Their numerical source
is [`../docs/COMAVI_current_results.md`](../docs/COMAVI_current_results.md).

The older top-level `COMAVI_Figure*.png` files, `manifest.json`, and scripts
under `src/` are retained as the historical pre-ISDS figure suite.


Publication figures for the COMAVI methods/benchmark manuscript, regenerated on the
canonical 61-variant benchmark (14 systems). See
`../docs/COMAVI_v7_canonical_benchmark_ledger.md` §12 for the recompute record and
`../reference_outputs/scored_61var_canonical.csv` for the underlying data.

| File | Manuscript figure | Content |
|------|-------------------|---------|
| `COMAVI_Figure1_pipeline_schematic.png`    | Figure 1 | Three-axis pipeline schematic (inputs → FoldX axes → concordance) |
| `COMAVI_Figure2_headline_competency.png`   | Figure 2 | Headline structural-agreement / mechanism-consistency threshold sweep + tier gradient |
| `COMAVI_Figure3_axis_competency.png`       | Figure 3 | Per-axis competency sweep + competency by mechanism class |
| `COMAVI_Figure4_measured_vs_foldx.png`     | Figure 4 | Predicted vs directly-measured ΔΔG (BRCA1-BRCT fold ρ=0.72; hemoglobin binding ρ=0.90) |
| `COMAVI_Figure5_alphamissense_vs_tier.png` | Figure 5 | AlphaMissense score vs COMAVI structural-evidence tier |

Headline (t = 2.5): mechanism-consistency 0.71, structural agreement 0.76 (92/121).
Tier gradient: Tier 1 100% → Tier 2 72% → Tier 3 70% → Tier 4 43% (Spearman ρ = −0.40, p = 0.0044).

## Regenerating

Each figure has a standalone script in [`src/`](src/). Run from the repository root:

```bash
python figures/src/figure1_pipeline_schematic.py
python figures/src/figure2_headline_competency.py
python figures/src/figure3_axis_competency.py
python figures/src/figure4_measured_vs_foldx.py
python figures/src/figure5_alphamissense_vs_tier.py
```

Scripts resolve all paths relative to the repository root and read only files
tracked here — no external inputs. Figures 2, 4 and 5 recompute their statistics
from `../reference_outputs/` and `../supplement/` at render time and print them, so
a run self-verifies against the frozen ledger values. Figures 1 and 3 render fixed
content (a schematic; the frozen sweep table) and print only the output path.
