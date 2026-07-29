#!/usr/bin/env python3
"""COMAVI Figure 3 — per-axis and per-mechanism-class competency across the
FoldX calling-threshold sweep (61-variant canonical benchmark).

Panel a: structural agreement per axis (monomer-fold / complex-fold / binding).
Panel b: mechanism-consistency per mechanism class.

Sweep numerators/denominators are the frozen canonical values recorded in
docs/COMAVI_v7_canonical_benchmark_ledger.md; this script renders them and does
not recompute the sweep (the sweep itself is run by the pipeline, see run.py).

Run from repo root:  python figures/src/figure3_axis_competency.py
"""
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "figures" / "COMAVI_Figure3_axis_competency.png"


def apply_figure_style(sizes=(9, 8, 7)):
    base, secondary, tick = sizes
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": base, "axes.labelsize": base,
        "axes.titlesize": base, "legend.fontsize": secondary,
        "xtick.labelsize": tick, "ytick.labelsize": tick, "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": 200, "savefig.dpi": 300,
        "savefig.bbox": "tight", "axes.titlelocation": "left",
        "lines.linewidth": 1.2, "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42})


def set_frame(ax):
    for s, v in zip(("top", "right", "bottom", "left"), (False, False, True, True)):
        ax.spines[s].set_visible(v)
    ax.tick_params(direction="out", length=3, width=0.6)


def panel_letter(ax, letter, dx=-0.18):
    ax.text(dx, 1.02, letter, transform=ax.transAxes, fontweight="bold",
            fontsize=10, va="bottom", ha="left")


TAGS = ["t1.0", "t1.5", "t2.0", "t2.5", "tSAP"]
XPOS = [1.0, 1.5, 2.0, 2.5, 3.2]
XLAB = ["1.0", "1.5", "2.0", "2.5", "SAP"]

C_MONO, C_CX, C_BIND = "#3a6ea5", "#2A7F62", "#b5495b"

# Panel a — direction-aware structural agreement, (correct, evaluable) per axis.
# Frozen canonical 61-variant sweep.
SWEEP = {
    "t1.0": {"monomer": (12, 16), "fold": (15, 25), "binding": (26, 32)},
    "t1.5": {"monomer": (13, 16), "fold": (20, 25), "binding": (26, 32)},
    "t2.0": {"monomer": (14, 16), "fold": (19, 25), "binding": (25, 32)},
    "t2.5": {"monomer": (14, 16), "fold": (20, 25), "binding": (24, 32)},
    "tSAP": {"monomer": (14, 16), "fold": (19, 25), "binding": (24, 32)},
}
AXINFO = [("monomer", "Monomer-fold (n=16)", C_MONO),
          ("fold", "Complex-fold (n=25)", C_CX),
          ("binding", "Binding / PPI (n=32)", C_BIND)]

# Panel b — mechanism-consistency by mechanism class, same threshold sweep.
CLASSB = {
    "Structurally-silent (n=31)": ("#7a7a7a", [0.419, 0.581, 0.677, 0.774, 0.806]),
    "Mixed-structural (n=9)":     ("#8e6fb0", [0.889, 0.889, 0.667, 0.667, 0.611]),
    "PPI-destabilization (n=6)":  ("#c26a2a", [0.667, 0.667, 0.667, 0.583, 0.583]),
    "Fold-mechanism \u00b7 BRCT (n=12)": ("#2a7f62", [0.667, 0.583, 0.750, 0.750, 0.667]),
}

apply_figure_style()
fig = plt.figure(figsize=(10.6, 4.7))
gs = fig.add_gridspec(1, 2, wspace=0.30, left=0.08, right=0.985, top=0.80, bottom=0.17)

ax1 = fig.add_subplot(gs[0, 0])
for key, lbl, col in AXINFO:
    vals = [SWEEP[t][key][0] / SWEEP[t][key][1] for t in TAGS]
    ax1.plot(XPOS[:4], vals[:4], "-o", color=col, lw=2.0, ms=5.5, zorder=3, label=lbl)
    ax1.plot([XPOS[3], XPOS[4]], [vals[3], vals[4]], "--", color=col, lw=1.2,
             alpha=0.8, zorder=2)
    ax1.plot(XPOS[4], vals[4], marker="*", color=col, ms=12, zorder=4,
             mec="white", mew=0.4)
ax1.axvline(2.5, color="#bbb", lw=1, ls=":", zorder=1)
ax1.set_xticks(XPOS); ax1.set_xticklabels(XLAB)
ax1.set_xlim(0.8, 3.55); ax1.set_ylim(0.55, 0.95)
ax1.set_xlabel("FoldX |\u0394\u0394G| call threshold (kcal/mol)")
ax1.set_ylabel("Structural agreement (direction-aware, per axis)")
ax1.set_title("Per-axis agreement across thresholds", loc="left")
ax1.text(2.55, 0.935, "canonical\nt=2.5", fontsize=6.5, color="#999", ha="left", va="top")
h, l = ax1.get_legend_handles_labels()
h += [Line2D([], [], marker="*", color="#666", ms=9, lw=0)]
l += ["SAP = Sapozhnikov CI"]
ax1.legend(h, l, loc="lower center", fontsize=6.6, handlelength=1.5,
           labelspacing=0.35, borderpad=0.5)
set_frame(ax1); panel_letter(ax1, "a")

ax2 = fig.add_subplot(gs[0, 1])
for lbl, (col, vals) in CLASSB.items():
    ax2.plot(XPOS[:4], vals[:4], "-o", color=col, lw=2.0, ms=5.0, zorder=3, label=lbl)
    ax2.plot([XPOS[3], XPOS[4]], [vals[3], vals[4]], "--", color=col, lw=1.2,
             alpha=0.8, zorder=2)
    ax2.plot(XPOS[4], vals[4], marker="*", color=col, ms=11, zorder=4,
             mec="white", mew=0.4)
ax2.axvline(2.5, color="#bbb", lw=1, ls=":", zorder=1)
ax2.set_xticks(XPOS); ax2.set_xticklabels(XLAB)
ax2.set_xlim(0.8, 3.55); ax2.set_ylim(0.30, 1.02)
ax2.set_xlabel("FoldX |\u0394\u0394G| call threshold (kcal/mol)")
ax2.set_ylabel("Mechanism-consistency")
ax2.set_title("Competency by mechanism class", loc="left")
ax2.text(2.55, 0.335, "canonical t=2.5", fontsize=6.5, color="#999", ha="left", va="bottom")
h2, l2 = ax2.get_legend_handles_labels()
h2 += [Line2D([], [], marker="*", color="#666", ms=9, lw=0)]
l2 += ["SAP = Sapozhnikov CI"]
ax2.legend(h2, l2, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=6.3,
           handlelength=1.4, labelspacing=0.3, borderpad=0.4)
set_frame(ax2); panel_letter(ax2, "b")

fig.suptitle("COMAVI benchmark competency (61-variant set; direction-aware "
             "agreement reproduces shipped-pipeline metric)",
             x=0.08, ha="left", fontsize=9.6, y=0.955)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
print("panel a t=2.5: monomer %d/%d, fold %d/%d, binding %d/%d"
      % (*SWEEP["t2.5"]["monomer"], *SWEEP["t2.5"]["fold"], *SWEEP["t2.5"]["binding"]))
print(f"saved {OUT.relative_to(REPO)}")
