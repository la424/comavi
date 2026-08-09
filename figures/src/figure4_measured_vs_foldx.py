#!/usr/bin/env python3
"""COMAVI Figure 4 — predicted vs. directly-measured DDG (two panels).

Panel a: fold axis, BRCA1 tandem-BRCT vs GdmCl unfolding (Rowling 2010).
Panel b: binding axis, hemoglobin tetramer W37 dose-series + N102T.

Data source: repository reference outputs (no hardcoded values).
Run from repo root:  python figures/src/figure4_measured_vs_foldx.py
"""
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "figures" / "COMAVI_Figure4_measured_vs_foldx.png"

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _comavi_style import (COLUMN_W as CW, BASE, SECOND, TICK, LETTER, panel_letter,
                           apply_figure_style, check_legibility, check_overlaps)

apply_figure_style()
mpl.rcParams.update({"font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"]})

C_FOLD_DESTAB = "#c0392b"
C_NEUTRAL = "#2980b9"
C_INTACT = "#8d8d8d"
C_BIND = "#b5495b"
C_FN = "#c0392b"
GREY = "#8a8a8a"

# ── Panel a data: BRCT fold axis ────────────────────────────────────────────
brct = pd.read_csv(REPO / "supplement/brct/brct_foldx_concordance.csv")
ROLE_STYLE = {
    "monomer_fold_destabilizer": ("Fold destabilizer", C_FOLD_DESTAB),
    "measured_neutral_control": ("Neutral control", C_NEUTRAL),
    "fold_intact_function_lost": ("Fold-intact / function-lost", C_INTACT),
}
# correlation is computed on the graded fold cohort only: the two
# fold-intact/function-lost variants (R1699L/Q) are excluded by design.
graded = brct[brct.role_in_cohort != "fold_intact_function_lost"]
rho_a, p_a = spearmanr(graded.measured_ddG_UF_kcal_mol, graded.foldx_ddg_monomer_mean)

# ── Panel b data: hemoglobin binding axis ───────────────────────────────────
conc = pd.read_csv(REPO / "reference_outputs/comavi_v7_concordance_annotated.csv",
                   low_memory=False)
HB_SERIES = ["W37Y", "W37A", "W37G", "W37E", "N102T"]
hbrows = conc[(conc.gene.astype(str).str.lower() == "hbb")
              & (conc.variant.isin(HB_SERIES))]
hb = pd.DataFrame({
    "variant": hbrows.variant.values,
    "measured": hbrows.expected_ddg_value.values,
    # alpha1-beta2 interface: the tetramer-assembly contact the measurements report
    "foldx": hbrows.ddg_binding_hba1_2.values,
})
rho_b, p_b = spearmanr(hb.measured, hb.foldx)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(CW, 3.5))

# ── Panel a ─────────────────────────────────────────────────────────────────
axA.plot([-2.5, 7.2], [-2.5, 7.2], ls=":", c="0.62", lw=0.9, zorder=1)
for role, (lab, col) in ROLE_STYLE.items():
    g = brct[brct.role_in_cohort == role]
    if g.empty:
        continue
    axA.errorbar(g.measured_ddG_UF_kcal_mol, g.foldx_ddg_monomer_mean,
                 yerr=g.foldx_ddg_sd, fmt="o", ms=7, mfc=col, mec="0.25",
                 mew=0.6, ecolor="0.55", elinewidth=0.8, capsize=2.2,
                 lw=0, label=lab, zorder=3)
# per-variant label offsets: the low-DDG cluster is crowded, so a few labels
# are placed away from the default upper-right position to avoid collisions.
LAB_OFF = {
    "L1664P": (1, -12), "R1751Q": (6, 4), "M1663K": (-6, 6),
    "R1699L": (6, -10), "V1665M": (6, -10), "P1806A": (6, 3),
    "A1843P": (6, 1), "Y1853C": (6, 1), "V1808A": (6, 2),
    "M1783T": (6, -2), "R1699Q": (6, 2),
}
LAB_HA = {"M1663K": "right"}
for _, r in brct.iterrows():
    if r.variant == "V1736A":
        continue
    axA.annotate(r.variant, (r.measured_ddG_UF_kcal_mol, r.foldx_ddg_monomer_mean),
                 textcoords="offset points", xytext=LAB_OFF.get(r.variant, (6, 2)),
                 ha=LAB_HA.get(r.variant, "left"),
                 fontsize=SECOND, color="#2b2b2b")

fn = brct[brct.variant == "V1736A"].iloc[0]
axA.annotate(
    f"FoldX false-negative\n(measured {fn.measured_ddG_UF_kcal_mol:.2f}, "
    f"FoldX {fn.foldx_ddg_monomer_mean:.2f})",
    xy=(fn.measured_ddG_UF_kcal_mol, fn.foldx_ddg_monomer_mean),
    xytext=(fn.measured_ddG_UF_kcal_mol + 0.6, fn.foldx_ddg_monomer_mean - 1.05),
    fontsize=SECOND, color=C_FN, ha="left", va="top",
    arrowprops=dict(arrowstyle="-", color=C_FN, lw=0.9))
axA.annotate("V1736A", (fn.measured_ddG_UF_kcal_mol, fn.foldx_ddg_monomer_mean),
             textcoords="offset points", xytext=(6, -10), fontsize=SECOND, color="#2b2b2b")

# Headroom for the rightmost/leftmost point labels, which are drawn in data
# space and would otherwise clip at the panel edge.
axA.set_xlim(-2.9, 8.1)
axA.set_xlabel("Measured $\\Delta\\Delta$G$_{U-F}$ (kcal/mol)")
axA.set_ylabel("FoldX $\\Delta\\Delta$G monomer (kcal/mol)")
axA.set_title(f"Fold axis — BRCA1 BRCT ({len(brct)} shown)", fontsize=BASE,
              loc="left", pad=4)
# Statistic sits inside the panel, not in the title: a two-line title at
# column width wraps into the axes.
axA.text(0.97, 0.03, f"Spearman \u03c1={rho_a:.2f}\nn={len(graded)}, p={p_a:.3f}",
         transform=axA.transAxes, ha="right", va="bottom",
         fontsize=SECOND, color="#333", linespacing=1.25)
axA.legend(loc="upper left", frameon=False, handletextpad=0.4,
           borderpad=0.2, labelspacing=0.3, fontsize=SECOND,
           bbox_to_anchor=(-0.01, 1.01))
panel_letter(axA, "a", dx=-0.20, dy=1.01)

# ── Panel b ─────────────────────────────────────────────────────────────────
axB.plot([-0.6, 10], [-0.6, 10], ls=":", c="0.62", lw=0.9, zorder=1)
axB.scatter(hb.measured, hb.foldx, s=52, c=C_BIND, edgecolors="0.25",
            linewidths=0.6, zorder=3, label="Hb tetramer assembly")
for _, r in hb.iterrows():
    axB.annotate(r.variant, (r.measured, r.foldx), textcoords="offset points",
                 xytext=(6, 2), fontsize=SECOND, color="#2b2b2b")
axB.set_xlim(-1.4, 12.2)
# x and y share the -2 gridline; suppress the redundant x tick at the corner
# so the two labels do not overprint.
axB.set_xticks([-2, 0, 2, 4, 6, 8, 10, 12][1:])
axB.set_xlabel("Measured assembly $\\Delta\\Delta$G (kcal/mol)")
axB.set_ylabel("FoldX $\\Delta\\Delta$G$_{bind}$ $\\alpha$1$\\beta$2 (kcal/mol)")
axB.set_title(f"Binding axis — hemoglobin tetramer (n={len(hb)})",
              fontsize=BASE, loc="left", pad=4)
axB.text(0.97, 0.03, f"Spearman \u03c1={rho_b:.2f}\np={p_b:.3f} \u00b7 W37 series",
         transform=axB.transAxes, ha="right", va="bottom",
         fontsize=SECOND, color="#333", linespacing=1.25)
axB.legend(loc="upper left", frameon=False, handletextpad=0.4, borderpad=0.2,
           fontsize=SECOND, bbox_to_anchor=(-0.01, 1.01))
panel_letter(axB, "b", dx=-0.20, dy=1.01)

fig.suptitle("Predicted vs directly-measured $\\Delta\\Delta$G",
             fontsize=BASE, y=1.0, x=0.012, ha="left", va="top",
             fontweight="bold")
# Source attributions move to the caption; at column width they do not fit
# under the axes without colliding with the panel below.
fig.tight_layout(rect=(0, 0, 1, 0.985), w_pad=1.5)
fig.subplots_adjust(top=0.845)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
_bad = check_legibility(fig)
_ov = check_overlaps(fig)
print(f"below-floor text: {_bad}")
print(f"overlaps: {_ov}")
print(f"panel a: rho={rho_a:.4f} p={p_a:.4f} n={len(graded)}")
print(f"panel b: rho={rho_b:.4f} p={p_b:.4f} n={len(hb)}")
print(f"saved {OUT.relative_to(REPO)}")
