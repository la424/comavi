"""Figure: threshold operating points.

Panels: (a) the recall/rejection trade as the FoldX calling threshold rises,
(b) how many experimentally measured destabilizers each threshold catches,
(c) which variants actually change grade when the threshold moves.

Provenance note
---------------
This figure was originally produced in an interactive session and existed
only in artifact lineage -- there was no script to re-run. The plotting code
was recovered from lineage and is committed here, re-laid-out at the shared
manuscript canvas width (see _comavi_style).
"""
import pathlib
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import apply_concordance_v5 as ac
from _comavi_style import (COLUMN_W as CW, BASE, SECOND, TICK,
                           apply_figure_style, set_frame, panel_letter,
                           check_legibility, check_overlaps)

OUT = REPO / "figures" / "COMAVI_Figure_operating_points.png"

GMAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
STRUCT = {"ppi_destab_mechanism", "mixed_structural", "fold_mechanism"}

# ---- data -----------------------------------------------------------------
df = pd.read_csv(REPO / "reference_outputs/scored_61var_canonical.csv")
partners = [p for p in ac.discover_partners(df) if f"ddg_{p}_confident" in df.columns]
tags = [t for t, _ in ac.THRESHOLD_SPECS]
unobs = set(ac.unobservable_variants())
g = df[~df.variant.isin(unobs)].copy()
for t in tags:
    g[f"g_{t}"] = g[f"mech_consistency_{t}"].map(GMAP)
com = g[g.expected_mech_class.isin(STRUCT | {"structurally_silent"})].copy()
com["truth_struct"] = com.expected_mech_class.isin(STRUCT)
com["max_abs_ddg"] = com.apply(lambda r: ac.compute_max_abs_ddg(r, partners), axis=1)

gc = [f"g_{t}" for t in tags]
recall = np.array([com[com.truth_struct][c].mean() for c in gc])
reject = np.array([com[~com.truth_struct][c].mean() for c in gc])
balanced = (recall + reject) / 2
n_struct = int(com.truth_struct.sum())
n_silent = int((~com.truth_struct).sum())

cal = pd.read_csv(REPO / "reference_outputs/COMAVI_delta_calibration_points.csv")
d = cal.dropna(subset=["measured_kcal", "foldx_ddg"]).copy()
d["abs_measured"] = d.measured_kcal.abs()
d["abs_foldx"] = d.foldx_ddg.abs()
REAL_DESTAB_KCAL = 1.0
real = d[d.abs_measured >= REAL_DESTAB_KCAL]
mf_rows = []
for tag, spec in ac.THRESHOLD_SPECS:
    tv = spec["fold"] if isinstance(spec, dict) else spec
    caught = real[real.abs_foldx >= tv]
    missed = real[real.abs_foldx < tv]
    mf_rows.append(dict(threshold=tag, foldx_t=tv, caught=len(caught),
                        missed=len(missed),
                        measured_recall=round(len(caught) / len(real), 3)))
mf = pd.DataFrame(mf_rows)
N_REAL = len(real)

com2 = com.copy()
com2["labile"] = com2[gc].nunique(axis=1) > 1
labdf = com2[["variant", "system", "truth_struct", "labile", "max_abs_ddg"]].copy()

# ---- figure ---------------------------------------------------------------
apply_figure_style()
x = np.arange(5)
XT = ["1.0", "1.5", "2.0", "2.5", "per-axis\nCI"]
C_R, C_S, C_G = "#C1440E", "#1F6FB4", "#555555"

fig, axes = plt.subplots(1, 3, figsize=(CW, 3.15))

ax = axes[0]
ax.plot(x, recall, "-o", color=C_R, lw=1.6, ms=4.0, zorder=4)
ax.plot(x, reject, "-s", color=C_S, lw=1.6, ms=4.0, zorder=4)
ax.plot(x, balanced, "--", color=C_G, lw=1.1, zorder=3)
ax.axvspan(1, 2, color="#F2C14E", alpha=0.22, zorder=0)
ax.text(1.5, 0.40, "crossover", color="#8A6D1E", fontsize=SECOND, ha="center", style="italic")
ax.text(-0.05, 0.80, f"recall on {n_struct}\nreal mechanisms", color=C_R,
        fontsize=SECOND, va="bottom", linespacing=1.2)
ax.text(4.15, 0.815, f"rejection of {n_silent}\nsilent variants", color=C_S,
        fontsize=SECOND, ha="right", va="top", linespacing=1.2)
ax.text(2.30, 0.655, "balanced", color=C_G, fontsize=SECOND, ha="left")
ax.annotate("canonical", xy=(3, 0.600), xytext=(2.05, 0.470), fontsize=SECOND,
            color="#333", arrowprops=dict(arrowstyle="->", lw=0.7, color="#333"))
ax.set_xticks(x); ax.set_xticklabels(XT)
ax.set_ylim(0.36, 0.95)
ax.set_xlabel("FoldX calling threshold (kcal/mol)")
ax.set_ylabel("Mechanism consistency")
ax.set_title("Raising the threshold trades\nrecall for rejection", loc="left", linespacing=1.25)
set_frame(ax)

ax = axes[1]
ax.bar(x, mf.caught, color=C_R, zorder=3, width=0.62)
ax.bar(x, mf.missed, bottom=mf.caught, color="#D9D9D9", zorder=3, width=0.62)
for i, r in mf.iterrows():
    ax.text(i, r.caught / 2, f"{int(r.caught)}", ha="center", va="center",
            fontsize=SECOND, color="w", fontweight="bold", zorder=4)
    ax.text(i, r.caught + r.missed / 2, f"{int(r.missed)}", ha="center",
            va="center", fontsize=SECOND, color="#555", zorder=4)
    ax.text(i, N_REAL + 2.0, f"{r.measured_recall:.2f}", ha="center",
            fontsize=SECOND, color=C_R)
ax.text(-0.60, N_REAL + 5.2, "recall on measured", fontsize=SECOND,
        color=C_R, ha="left", va="center")
ax.set_xticks(x); ax.set_xticklabels(XT)
ax.set_ylim(0, N_REAL + 8)
ax.set_xlabel("FoldX calling threshold (kcal/mol)")
ax.set_ylabel(f"Measured destabilizers\n(|ΔΔG| ≥ 1 kcal/mol, n = {N_REAL})")
ax.set_title(f"The canonical threshold misses\n{int(mf.missed.iloc[3])} of {N_REAL} confirmed destabilizers",
             loc="left", linespacing=1.25)
set_frame(ax)

ax = axes[2]
rngp = np.random.default_rng(3)
grps = [("stable", True, "Stable\nreal", C_R, 0.30),
        ("labile", True, "Labile\nreal", C_R, 1.0),
        ("labile", False, "Labile\nsilent", C_S, 1.0),
        ("stable", False, "Stable\nsilent", C_S, 0.30)]
YTOP = float(labdf.max_abs_ddg.max()) * 1.10
for i, (g_, ts, lab_, col, al) in enumerate(grps):
    s = labdf[(labdf.labile == (g_ == "labile")) & (labdf.truth_struct == ts)]
    ax.scatter(rngp.normal(i, 0.075, len(s)), s.max_abs_ddg, s=18, color=col,
               alpha=al, lw=0, zorder=3)
    m = s.max_abs_ddg.median()
    ax.plot([i - 0.24, i + 0.24], [m, m], color="k", lw=1.6, zorder=4)
    ax.text(i + 0.28, m, f"{m:.2f}", fontsize=SECOND, va="center")
    ax.text(i, YTOP - 0.5, f"n = {len(s)}", ha="center", va="top",
            fontsize=SECOND, color="#666")
ax.axhspan(1.0, 3.0, color="#F2C14E", alpha=0.18, zorder=0)
ax.annotate("threshold-\nsensitive band", xy=(2.0, 3.0), xytext=(2.0, 7.2),
            fontsize=SECOND, color="#8A6D1E", ha="center", va="bottom",
            style="italic", linespacing=1.2,
            arrowprops=dict(arrowstyle="-", lw=0.6, color="#8A6D1E"))
ax.set_xticks(range(4)); ax.set_xticklabels([q[2] for q in grps], linespacing=1.2)
ax.set_xlim(-0.55, 3.55)
ax.set_ylim(-0.4, YTOP)
ax.set_ylabel("Predicted max |ΔΔG| (kcal/mol)")
ax.set_title("Threshold choice only matters in a\nnarrow effect-size window",
             loc="left", linespacing=1.25)
set_frame(ax)

for a, l in zip(axes, "abc"):
    panel_letter(a, l)
fig.tight_layout(w_pad=1.6)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
print("below-floor text:", check_legibility(fig))
print("overlaps:", check_overlaps(fig))
