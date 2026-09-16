# Superseded figure files

The three figures below are tracked in this repository but are **not part of the
manuscript** and **cannot be regenerated**: no script in the repository produces
them. They were rendered before the v7.7 evidence-ledger correction, so the
numbers they display are stale.

| File | Status |
|---|---|
| `COMAVI_Figure5_scope.{png,pdf}` | no generator; superseded by `COMAVI_Figure5_alphamissense_vs_tier` |
| `COMAVI_Figure6_tier_null.{png,pdf}` | no generator; the manuscript's Figure 6 is `isds_v1/figure6_threshold_tradeoff.png` |
| `COMAVI_Figure_threshold_operating_points.{png,pdf}` | no generator; superseded by `COMAVI_Figure_operating_points` |

They are listed here rather than deleted so the decision is explicit and
reversible. To drop them from the release:

```
git rm figures/COMAVI_Figure5_scope.png figures/COMAVI_Figure5_scope.pdf \
       figures/COMAVI_Figure6_tier_null.png figures/COMAVI_Figure6_tier_null.pdf \
       figures/COMAVI_Figure_threshold_operating_points.png \
       figures/COMAVI_Figure_threshold_operating_points.pdf
```

Every other figure in `figures/` and `figures/isds_v1/` regenerates from a
committed generator; see `figures/isds_v1/SHA256SUMS.txt` for the checksummed
set and the repository README for the regeneration commands.
