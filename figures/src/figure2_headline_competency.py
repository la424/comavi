#!/usr/bin/env python3
"""COMAVI Figure 2 — headline competency and the structural-evidence tier gradient
(61-variant canonical benchmark).

Panel a: mechanism-consistency and structural agreement across the FoldX
  calling-threshold sweep (t = 1.0, 1.5, 2.0, 2.5, and the Sapozhnikov
  distinguishability criterion "Sap."), with the canonical operating point
  t = 2.5 marked.
Panel b: fraction pathogenic per structural-evidence tier, with Wilson score
  95% intervals.

Both panels recompute from reference_outputs/scored_61var_canonical.csv rather
than hard-coding, and print the computed values so a run self-verifies against
docs/COMAVI_v7_canonical_benchmark_ledger.md §17 (v7.3, current):
  MC/SA sweep  t1.0 0.579/0.718 · t1.5 0.649/0.756 · t2.0 0.684/0.756 ·
               t2.5 0.719/0.756 · tSAP 0.711/0.740   (SA denominator 131,
               MC graded n = 57 with the BRCA1-BRCT cohort pooled in)
  tier         T1 14/14 · T2 13/18 · T3 7/10 · T4 3/7;
               Fisher OR 3.78 p 0.080; Spearman rho -0.400 p 0.0044

Notes that are not evident from the code:
  * Structural agreement is aggregated as sum(numerator)/sum(denominator) over
    variants, NOT as the mean of per-variant ratios — the denominator counts
    gradeable AXES, so a per-variant mean would silently weight one-axis
    variants equally with four-axis ones.
  * The tier gradient is computed on the 49 tier-assigned rows only; the 12
    BRCT rows carry no tier (the tier is a PPI-interface construct) and drop out
    via dropna. This is why v7.3 pooling moved panel a but left panel b
    untouched — the cohort contributes mechanism and agreement rows, not tiers. "Pathogenic" includes pathogenic_gof — the tier measures
    structural-evidence strength, not direction of effect.
  * Fisher's exact test dichotomises tiers 1-2 vs 3-4; it is reported alongside
    the Spearman rho, which is the primary statistic (the dichotomy loses
    significance on the 61-set while the rank correlation survives).

Run from repo root:  python figures/src/figure2_headline_competency.py
"""
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

REPO = pathlib.Path(__file__).resolve().parents[2]
DATA = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "figures" / "COMAVI_Figure2_headline_competency.png"

META_GREY = "#888888"
C_MECH = "#b5495b"
C_SAGR = "#2A7F62"
C_TIER = "#3a6ea5"


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
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_letter(ax, letter, x=-0.13, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")


def wilson(k, n, z=1.96):
    """Wilson score interval — used instead of the normal approximation because
    Tier 1 is 14/14 and a Wald interval would collapse to zero width there."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


apply_figure_style()
canon = pd.read_csv(DATA, low_memory=False)

# ---- panel a: threshold sweep -------------------------------------------
MC_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
THRESH = ["t10", "t15", "t20", "t25", "tSAP"]
THRESH_LAB = ["1.0", "1.5", "2.0", "2.5", "Sap."]
CANON_I = 3                      # index of t = 2.5, the canonical operating point

mech, sagr = [], []
for t in THRESH:
    graded = canon[f"mech_consistency_{t}"].map(MC_MAP).dropna()
    num = int(canon[f"structural_agreement_n_{t}"].sum())
    den = int(canon[f"structural_agreement_d_{t}"].sum())
    mech.append(graded.mean())
    sagr.append(num / den)
    print(f"  {t:5s} MC={graded.mean():.3f} (n={len(graded)})  "
          f"SA={num/den:.3f} ({num}/{den})")

# ---- panel b: tier gradient ---------------------------------------------
tiered = canon.dropna(subset=["comavi_tier"]).copy()
tiered["tier_n"] = tiered["comavi_tier"].str.extract(r"(\d)").astype(int)
tiered["is_path"] = tiered["phenotype"].astype(str).str.startswith("pathogenic")

tiers = sorted(tiered["tier_n"].unique())
patho = [int(tiered.loc[tiered.tier_n == t, "is_path"].sum()) for t in tiers]
total = [int((tiered.tier_n == t).sum()) for t in tiers]
frac = np.array([k / n for k, n in zip(patho, total)])
lo = np.array([wilson(k, n)[0] for k, n in zip(patho, total)])
hi = np.array([wilson(k, n)[1] for k, n in zip(patho, total)])

strong = tiered.tier_n <= 2
table = [[int((strong & tiered.is_path).sum()), int((strong & ~tiered.is_path).sum())],
         [int((~strong & tiered.is_path).sum()), int((~strong & ~tiered.is_path).sum())]]
odds, p_fisher = fisher_exact(table)
rho, p_rho = spearmanr(tiered.tier_n, tiered.is_path.astype(int))
print(f"  tier {list(zip(patho, total))}  OR={odds:.2f} p={p_fisher:.3f}  "
      f"rho={rho:.3f} p={p_rho:.4f}  (n={len(tiered)})")

# ---- render --------------------------------------------------------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.1))

x = np.arange(len(THRESH_LAB))
axA.plot(x, mech, "-o", color=C_MECH, lw=1.8, ms=5, zorder=3)
axA.plot(x, sagr, "-s", color=C_SAGR, lw=1.8, ms=4.5, zorder=3)
axA.axvline(CANON_I, color=META_GREY, lw=0.8, ls=(0, (3, 3)), zorder=1)
axA.set_xticks(x)
axA.set_xticklabels(THRESH_LAB)
axA.set_xlabel("FoldX ΔΔG threshold (kcal/mol)")
axA.set_ylabel("Fraction correct")
axA.set_ylim(0.45, 0.85)
axA.margins(x=0.05)
# series labelled at the right-hand end rather than in a legend box, so the
# reader does not have to map marker shapes back to names
axA.text(x[-1] + 0.08, sagr[-1] + 0.006, "Structural\nagreement",
         color=C_SAGR, va="center", ha="left", fontsize=7, linespacing=0.95)
axA.text(x[-1] + 0.08, mech[-1] - 0.012, "Mechanism-\nconsistency",
         color=C_MECH, va="center", ha="left", fontsize=7, linespacing=0.95)
axA.annotate("canonical\nt = 2.5", xy=(CANON_I, 0.47), ha="center", va="bottom",
             fontsize=6.5, color=META_GREY)
axA.set_title("Mechanism calls match literature ground truth", fontsize=8)
axA.set_xlim(-0.2, len(x) - 0.35 + 1.15)   # right margin holds the end labels
set_frame(axA)
panel_letter(axA, "a")

xb = np.arange(len(tiers))
axB.bar(xb, frac, width=0.66, color=C_TIER, edgecolor="white", lw=0.5, zorder=2)
axB.errorbar(xb, frac, yerr=np.vstack([frac - lo, hi - frac]), fmt="none",
             ecolor="#33333399", elinewidth=1.0, capsize=2.5, zorder=3)
for xi, (f, k, n) in enumerate(zip(frac, patho, total)):
    axB.text(xi, min(f + (hi - frac)[xi] + 0.03, 1.02), f"{k}/{n}",
             ha="center", va="bottom", fontsize=6.8)
axB.set_xticks(xb)
axB.set_xticklabels([str(t) for t in tiers])
axB.set_xlabel("Structural-evidence tier")
axB.set_ylabel("Fraction pathogenic")
axB.set_ylim(0, 1.14)
axB.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
for xi, lab in ((0, "strongest"), (len(tiers) - 1, "weakest")):
    axB.annotate(lab, xy=(xi, 0.0), xytext=(xi, -0.20),
                 textcoords=axB.get_xaxis_transform(), ha="center", va="top",
                 fontsize=6.3, color=META_GREY)
axB.text(0.5, 0.965,
         f"OR = {odds:.1f}, p = {p_fisher:.2f}; ρ = {rho:.2f}, p = {p_rho:.3f}"
         .replace("-", "\u2212"),
         transform=axB.transAxes, ha="center", va="top", fontsize=6.0,
         color="#222")
axB.set_title("Tier is a monotonic pathogenicity gradient", fontsize=8)
set_frame(axB)
panel_letter(axB, "b")

fig.subplots_adjust(left=0.085, right=0.985, top=0.86, bottom=0.20, wspace=0.34)
fig.savefig(OUT, dpi=300)
fig.savefig(OUT.with_suffix(".pdf"))
print(f"saved {OUT.relative_to(REPO)}")
