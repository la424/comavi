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

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "savefig.dpi": 300, "figure.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
})

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

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.6, 4.9))

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
    "L1664P": (2, -11), "R1751Q": (7, 4), "M1663K": (7, -3),
    "R1699L": (-40, -3), "V1665M": (7, -9), "P1806A": (7, 2),
}
for _, r in brct.iterrows():
    if r.variant == "V1736A":
        continue
    axA.annotate(r.variant, (r.measured_ddG_UF_kcal_mol, r.foldx_ddg_monomer_mean),
                 textcoords="offset points", xytext=LAB_OFF.get(r.variant, (7, 2)),
                 fontsize=7, color="#2b2b2b")

fn = brct[brct.variant == "V1736A"].iloc[0]
axA.annotate(
    f"FoldX false-negative\n(measured {fn.measured_ddG_UF_kcal_mol:.2f}, "
    f"FoldX {fn.foldx_ddg_monomer_mean:.2f})",
    xy=(fn.measured_ddG_UF_kcal_mol, fn.foldx_ddg_monomer_mean),
    xytext=(fn.measured_ddG_UF_kcal_mol + 0.6, fn.foldx_ddg_monomer_mean - 1.05),
    fontsize=7, color=C_FN, ha="left", va="top",
    arrowprops=dict(arrowstyle="-", color=C_FN, lw=0.9))
axA.annotate("V1736A", (fn.measured_ddG_UF_kcal_mol, fn.foldx_ddg_monomer_mean),
             textcoords="offset points", xytext=(7, -9), fontsize=7, color="#2b2b2b")

axA.set_xlabel(r"Measured $\Delta\Delta$G$_{U-F}$ (kcal/mol) $\cdot$ Rowling 2010 GdmCl unfolding")
axA.set_ylabel(r"COMAVI FoldX $\Delta\Delta$G monomer (kcal/mol)")
axA.set_title(f"Fold axis — BRCA1 BRCT ({len(brct)} shown)\n"
              f"Spearman \u03c1={rho_a:.2f} (n={len(graded)}, p={p_a:.3f}); "
              f"R1699L/Q excluded (fold-intact)", fontsize=8.6, loc="left")
axA.legend(loc="upper left", frameon=False, handletextpad=0.4,
           borderpad=0.2, labelspacing=0.35, bbox_to_anchor=(0.0, 0.99))
axA.text(-0.115, 1.06, "a", transform=axA.transAxes, fontsize=13,
         fontweight="bold", va="top")

# ── Panel b ─────────────────────────────────────────────────────────────────
axB.plot([-0.6, 10], [-0.6, 10], ls=":", c="0.62", lw=0.9, zorder=1)
axB.scatter(hb.measured, hb.foldx, s=52, c=C_BIND, edgecolors="0.25",
            linewidths=0.6, zorder=3, label="Hb tetramer assembly")
for _, r in hb.iterrows():
    axB.annotate(r.variant, (r.measured, r.foldx), textcoords="offset points",
                 xytext=(7, 2), fontsize=7, color="#2b2b2b")
axB.set_xlabel(r"Measured assembly $\Delta\Delta$G (kcal/mol) $\cdot$ tetramer$\rightarrow$dimer"
               "\nKwiatkowski 1998 / Bonaventura 1968")
axB.set_ylabel(r"COMAVI FoldX $\Delta\Delta$G$_{bind}$ $\alpha$1$\beta$2 (kcal/mol)")
axB.set_title(f"Binding axis — hemoglobin tetramer (n={len(hb)})\n"
              f"Spearman \u03c1={rho_b:.2f} (p={p_b:.3f}) \u00b7 W37 dose-series",
              fontsize=8.6, loc="left")
axB.legend(loc="upper left", frameon=False, handletextpad=0.4, borderpad=0.2)
axB.text(-0.115, 1.06, "b", transform=axB.transAxes, fontsize=13,
         fontweight="bold", va="top")

fig.suptitle("COMAVI predicted vs directly-measured "
             r"$\Delta\Delta$G — the two systems with quantitative biophysical comparators",
             fontsize=9.6, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.965), w_pad=2.6)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"panel a: rho={rho_a:.4f} p={p_a:.4f} n={len(graded)}")
print(f"panel b: rho={rho_b:.4f} p={p_b:.4f} n={len(hb)}")
print(f"saved {OUT.relative_to(REPO)}")
