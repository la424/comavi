#!/usr/bin/env python3
"""COMAVI Figure S2 — binding axis vs. directly-measured DDG on two SKEMPI complexes.

Panel a: barnase-barstar (PDB 1BRS), 24 substitutions. The seven buried
         active-site barnase Glu73 substitutions are measured destabilizing but
         predicted stabilizing (FoldX electrostatics failure mode) and are drawn
         as a visually distinct series, excluded from the reported rho.
Panel b: TEM1 beta-lactamase-BLIP (PDB 1JTG), BLIP side only, 24 substitutions.

Cross-system corroboration of the within-system hemoglobin anchor (Figure 4b).
None of these values enters any main-text benchmark metric.

Data source: repository supplement tables (no hardcoded statistics).
Run from repo root:  python figures/src/figureS2_skempi_validation.py
"""
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr, pearsonr

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "figures" / "COMAVI_FigureS2_skempi_validation.png"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "savefig.dpi": 300, "figure.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
})

C_IFACE = "#2980b9"     # interface substitutions (same blue as neutral-control in Fig 4)
C_EXCL = "#c0392b"      # excluded Glu73 cluster (same red as fold-destabilizer in Fig 4)
GREY = "#8a8a8a"

bb = pd.read_csv(REPO / "supplement/skempi/bb_binding_validation.csv")
jt = pd.read_csv(REPO / "supplement/skempi/jt_binding_validation.csv")

# Glu73 cluster: FoldX code E<chain>73<mut> on the barnase chain.
g73 = bb["foldx_code"].str.match(r"E[A-Z]?73")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))

def identity_line(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], ls="--", lw=0.9, color=GREY, zorder=1)

def frame(ax, x, y, pad=0.10):
    lo = min(x.min(), y.min()); hi = max(x.max(), y.max())
    span = hi - lo
    ax.set_xlim(lo - pad * span, hi + pad * span)
    ax.set_ylim(lo - pad * span, hi + pad * span)
    return lo - pad * span, hi + pad * span

# ── Panel a: barnase-barstar ────────────────────────────────────────────────
ax = axes[0]
keep, drop = bb[~g73], bb[g73]
lo, hi = frame(ax, bb["ddg_meas"], bb["ddg_pred_bind"])
identity_line(ax, lo, hi)
ax.axhline(0, lw=0.6, color="#d8d8d8", zorder=0)
ax.axvline(0, lw=0.6, color="#d8d8d8", zorder=0)
ax.scatter(keep["ddg_meas"], keep["ddg_pred_bind"], s=42, facecolor=C_IFACE,
           edgecolor="white", linewidth=0.6, zorder=3,
           label=f"interface (n={len(keep)})")
ax.scatter(drop["ddg_meas"], drop["ddg_pred_bind"], s=42, marker="s",
           facecolor="none", edgecolor=C_EXCL, linewidth=1.2, zorder=3,
           label=f"buried Glu73, excluded (n={len(drop)})")
rho_a, p_a = spearmanr(keep["ddg_meas"], keep["ddg_pred_bind"])
ax.set_title("Barnase\u2013barstar (1BRS)")
ax.set_xlabel("Measured $\\Delta\\Delta G_{bind}$ (kcal/mol)")
ax.set_ylabel("FoldX predicted $\\Delta\\Delta G_{bind}$ (kcal/mol)")
ax.annotate(f"$\\rho$ = {rho_a:.2f} (p = {p_a:.3f})\nexcluding Glu73 cluster",
            xy=(0.03, 0.97), xycoords="axes fraction", va="top", ha="left",
            fontsize=7.5, linespacing=1.35)
ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.83), frameon=False,
          handletextpad=0.4, labelspacing=0.3)

# ── Panel b: TEM1-BLIP ─────────────────────────────────────────────────────
ax = axes[1]
lo, hi = frame(ax, jt["ddg_meas"], jt["ddg_pred_bind"])
identity_line(ax, lo, hi)
ax.axhline(0, lw=0.6, color="#d8d8d8", zorder=0)
ax.axvline(0, lw=0.6, color="#d8d8d8", zorder=0)
ax.scatter(jt["ddg_meas"], jt["ddg_pred_bind"], s=42, facecolor=C_IFACE,
           edgecolor="white", linewidth=0.6, zorder=3)
rho_b, p_b = spearmanr(jt["ddg_meas"], jt["ddg_pred_bind"])
ax.set_title("TEM1\u2013BLIP, BLIP side (1JTG)")
ax.set_xlabel("Measured $\\Delta\\Delta G_{bind}$ (kcal/mol)")
ax.set_ylabel("FoldX predicted $\\Delta\\Delta G_{bind}$ (kcal/mol)")
ax.annotate(f"$\\rho$ = {rho_b:.2f} (p = {p_b:.3f})\n"
            f"interface (n={len(jt)}), no exclusions",
            xy=(0.03, 0.97), xycoords="axes fraction", va="top", ha="left",
            fontsize=7.5, linespacing=1.35)

for letter, ax in zip("ab", axes):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="left")

fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")

# ── Reported statistics (stdout: the numbers the caption and Note S1 cite) ──
def stats(df, lab):
    r, pr = pearsonr(df["ddg_meas"], df["ddg_pred_bind"])
    rho, prho = spearmanr(df["ddg_meas"], df["ddg_pred_bind"])
    d = df[df["ddg_meas"] > 0]
    nsc = int((d["ddg_pred_bind"] > 0).sum())
    print(f"  {lab:<34s} n={len(df):2d}  r={r:.2f} (p={pr:.4f})  "
          f"rho={rho:.2f} (p={prho:.4f})  sign-correct={nsc}/{len(d)} "
          f"({100 * nsc / len(d):.0f}%)")

print(f"wrote {OUT.name} and {OUT.with_suffix('.pdf').name}")
stats(bb, "Barnase-barstar, all")
stats(bb[~g73], "Barnase-barstar, ex-Glu73")
stats(jt, "TEM1-BLIP (BLIP side)")
comb = pd.concat([bb[~g73][["ddg_meas", "ddg_pred_bind"]],
                  jt[["ddg_meas", "ddg_pred_bind"]]])
stats(comb, "Combined, ex-Glu73")
