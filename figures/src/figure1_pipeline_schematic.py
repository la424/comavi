#!/usr/bin/env python3
"""COMAVI Figure 1 — pipeline schematic.

Pure schematic: no data inputs, no computed values. Layout is in a 0-100
coordinate box so element positions read as percentages of the canvas.

Design notes worth preserving:
  * The "both inputs" italic label sits under the upper input box, to the LEFT
    of the elbow connector. It must not sit on the vertical merge line.
  * The three FoldX axes are drawn inside one group box, because the figure's
    point is that they are three readouts of one structural model, not three
    pipelines.
  * The structural-evidence tier is a separate lane feeding concordance from
    below, to show it is computed without any DDG term (non-circular).

Run from repo root:  python figures/src/figure1_pipeline_schematic.py
"""
import pathlib

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "figures" / "COMAVI_Figure1_pipeline_schematic.png"

# DejaVu Sans first (not Helvetica): the schematic uses circled digits and an
# arrow glyph that Helvetica lacks, and font fallback renders them unevenly.
plt.rcParams.update({"font.size": 8, "font.family": "sans-serif",
                     "font.sans-serif": ["DejaVu Sans"]})

C_IN = "#5b6b7a"      # inputs, slate
C_MONO = "#3a6ea5"    # monomer-fold axis, blue
C_CX = "#2A7F62"      # complex-fold axis, green
C_BIND = "#b5495b"    # binding axis, rose
C_TIER = "#c8912b"    # tier lane, amber
C_CONC = "#3d4451"    # concordance, dark slate
C_OUT = "#2f2f2f"
C_AMB_T = "#8a6410"   # amber text
C_AMB_B = "#5a4708"   # amber body text


def tint(hexc, a=0.14):
    """Very light fill derived from a stroke colour."""
    r = int(hexc[1:3], 16) / 255
    g = int(hexc[3:5], 16) / 255
    b = int(hexc[5:7], 16) / 255
    return (1 - a * (1 - r), 1 - a * (1 - g), 1 - a * (1 - b))


fig, ax = plt.subplots(figsize=(13.8, 7.8))
ax.set_xlim(0, 100)
ax.set_ylim(21, 95)   # crop empty canvas below the tier lane
ax.axis("off")


def box(x, y, w, h, title, body="", ec=C_IN, tfs=9.2, bfs=7.6, fc=None,
        tcol=None, lw=1.6, ls="-", body_dy=7.0, bcol="#1a1a1a", bstyle=None):
    fc = tint(ec) if fc is None else fc
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.35,rounding_size=1.6",
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=2))
    cx = x + w / 2
    if body:
        ax.text(cx, y + h - 2.4, title, ha="center", va="top", fontsize=tfs,
                fontweight="bold", color=tcol or ec, zorder=3)
        ax.text(cx, y + h - body_dy, body, ha="center", va="top", fontsize=bfs,
                color=bcol, zorder=3, linespacing=1.35, style=bstyle)
    else:
        ax.text(cx, y + h / 2, title, ha="center", va="center", fontsize=tfs,
                fontweight="bold", color=tcol or ec, zorder=3, linespacing=1.3)


def arrow(x1, y1, x2, y2, col="#555", lw=1.8, rad=0.0, ls="-", style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=15, lw=lw, color=col, zorder=1,
                                 connectionstyle=f"arc3,rad={rad}", linestyle=ls))


# ---------------- STAGE A: inputs ----------------
IN_X, IN_W = 2.3, 16.6
box(IN_X, 59.8, IN_W, 13.1, "Missense variant",
    "protein \u00b7 residue\nWT \u2192 Mut", ec=C_IN, body_dy=7.4)
box(IN_X, 40.3, IN_W, 12.7, "Structure",
    "AlphaFold-3 multimer\nor experimental PDB", ec=C_IN, body_dy=7.4)
ax.text(IN_X + IN_W / 2, 75.8, "INPUTS", ha="center", fontsize=8.5,
        fontweight="bold", color=C_IN, alpha=0.8)

# "both inputs": placed under the UPPER input box, clear of the merge line at
# ELBOW_X. Do not centre it on the elbow -- it collides with the vertical line.
ax.text(IN_X + IN_W / 2, 56.3, "both inputs", ha="center", va="center",
        fontsize=6.8, style="italic", color="#7a8794", zorder=4)

# elbow: both inputs merge onto one vertical line, then one junction point
ELBOW_X, JUNC_Y = 21.6, 56.4
Y_TOP, Y_BOT = 66.3, 46.6
for yy in (Y_TOP, Y_BOT):
    ax.plot([IN_X + IN_W, ELBOW_X], [yy, yy], color=C_IN, lw=1.4, zorder=1)
ax.plot([ELBOW_X, ELBOW_X], [Y_BOT, Y_TOP], color=C_IN, lw=1.4, zorder=1)
ax.plot([ELBOW_X], [JUNC_Y], marker="o", ms=5, color=C_IN, zorder=3)

# ---------------- STAGE B: three FoldX axes ----------------
GX, GY, GW, GH = 25.5, 48.9, 27.8, 42.5
ax.add_patch(FancyBboxPatch((GX, GY), GW, GH,
                            boxstyle="round,pad=0.4,rounding_size=1.6",
                            fc="#fafbfc", ec="#c6ced6", lw=1.2, zorder=1))
gtitle = "Three independent biophysical axes  (FoldX 5.1)"
ax.text(GX + GW / 2, GY + GH - 2.0, gtitle, ha="center", va="top",
        fontsize=9.4, fontweight="bold", color="#1a1a1a", zorder=3)
# underline the group title
ax.plot([GX + 2.6, GX + GW - 2.6], [GY + GH - 3.9] * 2, color="#1a1a1a",
        lw=0.9, zorder=3)

AX_X, AX_W = 26.6, 25.6
box(AX_X, 78.5, AX_W, 8.5, "\u2460 Monomer-fold stability",
    "BuildModel \u0394\u0394G \u00b7 isolated subunit", ec=C_MONO, tfs=9.0,
    bfs=7.5, body_dy=6.4)
box(AX_X, 68.6, AX_W, 8.4, "\u2461 Complex-fold stability",
    "BuildModel \u0394\u0394G \u00b7 within complex", ec=C_CX, tfs=9.0,
    bfs=7.5, body_dy=6.4)
box(AX_X, 58.5, AX_W, 8.5, "\u2462 Binding-interface disruption",
    "AnalyseComplex \u0394\u0394G$_{bind}$ \u00b7 top partner", ec=C_BIND,
    tfs=9.0, bfs=7.5, body_dy=6.4)

# pLDDT gate: one compact dashed band inside the group box
GTX, GTY, GTW, GTH = AX_X, 51.1, AX_W, 4.2
ax.add_patch(FancyBboxPatch((GTX, GTY), GTW, GTH,
                            boxstyle="round,pad=0.3,rounding_size=1.2",
                            fc="#fdf7ea", ec=C_TIER, lw=1.2, ls="--", zorder=2))
ax.text(GTX + GTW / 2, GTY + GTH / 2,
        "pLDDT gate \u00b7 contacts \u226550 \u00b7 \u00d70.7 (50\u201370), "
        "\u00d70.4 (<50) \u00b7 exp. bypass",
        ha="center", va="center", fontsize=6.5, fontweight="bold",
        color=C_AMB_T, zorder=3)

# ---------------- STAGE C: independent tier lane ----------------
TL_X, TL_Y, TL_W, TL_H = 25.5, 25.0, 27.8, 16.0
ax.add_patch(FancyBboxPatch((TL_X, TL_Y), TL_W, TL_H,
                            boxstyle="round,pad=0.4,rounding_size=1.6",
                            fc=tint(C_TIER), ec=C_TIER, lw=1.9, zorder=2))
ax.text(TL_X + TL_W / 2, TL_Y + TL_H - 2.1,
        "Structural-evidence tier  (Tier 1 \u2192 4)", ha="center", va="top",
        fontsize=9.2, fontweight="bold", color=C_AMB_T, zorder=3)
ax.text(TL_X + TL_W / 2, TL_Y + TL_H - 5.6,
        "independent, non-circular yardstick", ha="center", va="top",
        fontsize=7.4, style="italic", color=C_AMB_T, zorder=3)
ax.text(TL_X + TL_W / 2, TL_Y + TL_H - 8.6,
        "composite disruption score = Grantham severity \u00d7\n"
        "contacts (intra-chain + multimer interface) + burial bonus\n"
        "+ pLDDT multiplier   \u00b7   NO \u0394\u0394G term",
        ha="center", va="top", fontsize=7.0, color=C_AMB_B,
        linespacing=1.45, zorder=3)

# junction -> axes group, and junction -> tier lane
arrow(ELBOW_X + 0.3, JUNC_Y, GX - 0.6, 69.5, col=C_IN, lw=1.6, rad=-0.16)
arrow(ELBOW_X + 0.3, JUNC_Y, TL_X - 0.6, 35.0, col=C_IN, lw=1.6, rad=0.16)

# ---------------- STAGE D: concordance ----------------
CC_X, CC_Y, CC_W, CC_H = 57.9, 44.8, 22.5, 35.3
ax.add_patch(FancyBboxPatch((CC_X, CC_Y), CC_W, CC_H,
                            boxstyle="round,pad=0.35,rounding_size=1.6",
                            fc=tint(C_CONC), ec=C_CONC, lw=2.0, zorder=2))
# title pinned to the top of the box (not vertically centred) so it clears the
# four-item comparator list below it
ax.text(CC_X + CC_W / 2, CC_Y + CC_H - 2.6, "Four-way concordance",
        ha="center", va="top", fontsize=9.6, fontweight="bold", color=C_CONC,
        zorder=3)
ax.text(CC_X + CC_W / 2, CC_Y + CC_H - 8.0, "COMAVI structural call  vs.",
        ha="center", va="top", fontsize=8.0, style="italic", color="#1a1a1a",
        zorder=3)
items = [("Structural-evidence tier", C_TIER),
         ("FoldX \u0394\u0394G", C_CONC),
         ("AlphaMissense", C_CONC),
         ("Clinical (Franklin)", C_CONC)]
for i, (t, c) in enumerate(items):
    yy = CC_Y + CC_H - 13.0 - i * 4.5
    ax.text(CC_X + 2.6, yy, "\u2022", ha="left", va="center", fontsize=11,
            color=c, fontweight="bold", zorder=3)
    ax.text(CC_X + 5.2, yy, t, ha="left", va="center", fontsize=7.8,
            color="#1a1a1a", zorder=3)

arrow(GX + GW + 0.3, 73.0, CC_X - 0.6, CC_Y + CC_H - 5.0, col="#6b7783",
      lw=1.7, rad=-0.06)
arrow(TL_X + TL_W + 0.3, TL_Y + TL_H / 2, CC_X + 2.5, CC_Y - 0.5, col=C_TIER,
      lw=1.9, rad=-0.24)

# ---------------- STAGE E: outputs ----------------
box(84.0, 59.8, 15.4, 17.6, "MECHANISM",
    "which structural axis\nis disrupted:\nfold, complex,\nor binding interface",
    ec=C_BIND, tfs=9.4, bfs=7.4, tcol=C_BIND, body_dy=7.6)
box(84.0, 44.2, 15.4, 13.7, "STRENGTH",
    "how severe the\ndisruption is:\nFoldX \u0394\u0394G magnitude,\ngraded by threshold",
    ec=C_MONO, tfs=9.4, bfs=7.4, tcol=C_MONO, body_dy=6.4)
ax.text(91.7, 80.0, "OUTPUT", ha="center", fontsize=8.5, fontweight="bold",
        color=C_OUT, alpha=0.8)

arrow(CC_X + CC_W + 0.3, CC_Y + CC_H - 6.0, 83.6, 70.0, col=C_CONC, lw=1.8,
      rad=0.04)
arrow(CC_X + CC_W + 0.3, CC_Y + 9.0, 83.6, 52.5, col=C_CONC, lw=1.8, rad=-0.04)

# the deliberate separation of mechanism/strength from any clinical verdict
ax.add_patch(FancyBboxPatch((82.5, 35.7), 18.4, 5.5,
                            boxstyle="round,pad=0.3,rounding_size=1.2",
                            fc="white", ec="#999", lw=1.0, ls=":", zorder=2))
ax.text(91.7, 38.45, "reported separately from\nany pathogenicity verdict",
        ha="center", va="center", fontsize=7.2, style="italic", color="#444",
        linespacing=1.25, zorder=3)

fig.tight_layout(pad=0.3)
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
print(f"saved {OUT.relative_to(REPO)}")
