# COMAVI main figures (methods/benchmark paper, 61-variant set)

<!-- COMAVI_ISDS_V1_FIGURE_BANNER -->

## The manuscript's figures are under `submission/`

[`submission/`](submission/) holds the figures the manuscript embeds and the
ones uploaded with it, built by
[`submission/make_submission_figures.py`](submission/make_submission_figures.py)
and gated by `scripts/verify_embedded_figures.py`,
`scripts/verify_figure_specs.py` and `scripts/verify_si_package.py`. Figure
numbering follows the manuscript's own `@si` caption block, which is the
authority.

Everything else in this directory is **superseded and carries pre-correction
numbers**:

| Location | Status |
|---|---|
| `submission/` | current — the manuscript's figures |
| `isds_v1/` | superseded ISDS-v1 development suite |
| top-level `COMAVI_Figure*.png`, `manifest.json`, `src/` | superseded pre-ISDS suite |

The superseded files are retained as a development record. **They are not the
manuscript's figures and are deliberately not mapped to its figure numbers**,
because that mapping was the defect this section replaces: a table here
assigned "Figure 1" to "Figure 5" to the top-level PNGs and quoted headline
values beside them, and both went stale. The quoted values were a third copy of
numbers the manuscript and the root `README.md` already carry, so they are
removed rather than updated — a duplicated number with no generator is what
produced the drift. For current values read `../reference_outputs/` or the
Results section; both are regenerated and gated.

See `../docs/COMAVI_v7_canonical_benchmark_ledger.md` for the recompute record
and `../reference_outputs/scored_61var_canonical.csv` for the underlying data.

## Regenerating

For the manuscript's figures, run from the repository root:

```bash
python figures/submission/make_submission_figures.py
python scripts/embed_figures.py --docx docs/COMAVI_manuscript_current.docx
python scripts/verify_embedded_figures.py --docx docs/COMAVI_manuscript_current.docx
python scripts/verify_figure_specs.py --render
python scripts/verify_si_package.py
```

The commands below rebuild the **superseded** pre-ISDS suite. They are kept so
the development record stays reproducible; their output is not the manuscript's
figures and will not match its numbers.

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
