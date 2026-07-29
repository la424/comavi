#!/usr/bin/env python3
"""COMAVI Figure 5 — COMAVI structural-evidence tier vs. AlphaMissense score.

Shows that the structural-evidence tier and a sequence-based pathogenicity
predictor carry partially independent signal: pathogenic variants appear across
all four tiers, and tier-1 variants span the full AlphaMissense range.

Inputs (repository reference outputs):
  reference_outputs/scored_61var_canonical.csv
  supplement/brct/scored_56var_with_brct.csv   (phenotype labels)

Run from repo root:  python figures/src/figure5_alphamissense_vs_tier.py
"""
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "figures" / "COMAVI_Figure5_alphamissense_vs_tier.png"


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
        "lines.linewidth": 1.2, "pdf.fonttype": 42, "ps.fonttype": 42})


apply_figure_style()

canon = pd.read_csv(REPO / "reference_outputs/scored_61var_canonical.csv", low_memory=False)
bv = pd.read_csv(REPO / "supplement/brct/scored_56var_with_brct.csv", low_memory=False)
bv["key"] = list(zip(bv["gene"].str.lower(), bv["variant"]))
pmap = dict(zip(bv["key"], bv["phenotype"]))

# phenotype labels for the two systems added after the 56-variant set
new_phen = {("cfh", "R78G"): "pathogenic", ("cfh", "R53H"): "pathogenic",
            ("cfh", "I62V"): "benign",
            ("vwf", "R1334Q"): "pathogenic", ("vwf", "A1381T"): "benign"}


def role(k):
    return new_phen.get(k, pmap.get(k))


# BRCA1-BRCT is a monomer-fold system with no PPI axis; excluded from this panel
ppi = canon[canon["system"] != "brca1_brct"].copy()
ppi["key"] = list(zip(ppi["gene"].str.lower(), ppi["variant"]))
ppi["role"] = ppi["key"].map(role)
ppi["_ispath"] = ppi["role"].astype(str).str.contains("patho")
ppi["_tier"] = ppi["comavi_tier"].str.extract(r"(\d)").astype(int)
ppi["AM_score"] = pd.to_numeric(ppi["AM pathogenicity"], errors="coerce")

pv = ppi[ppi["AM_score"].notna()].copy()
npath = int(pv["_ispath"].sum())
nben = len(pv) - npath

C_PATH, C_BEN = "#b5495b", "#3a6ea5"
rng = np.random.default_rng(7)          # fixed seed: jitter is reproducible
pv["_y"] = pv["_tier"] + rng.uniform(-0.16, 0.16, len(pv))

fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.axvspan(0.34, 0.564, color="0.88", zorder=0)
ax.text(0.452, 0.60, "AM\nambiguous", ha="center", va="center", fontsize=7,
        color="0.45", style="italic")
for isp, col, lab in [(True, C_PATH, f"Pathogenic (n={npath})"),
                      (False, C_BEN, f"Benign (n={nben})")]:
    s = pv[pv["_ispath"] == isp]
    ax.scatter(s["AM_score"], s["_y"], c=col, s=46, alpha=0.85,
               edgecolor="white", linewidth=0.6, zorder=3, label=lab)

div = {("pik3ca", "H1047R"): ("PIK3CA H1047R", 10, 14),
       ("tnni3", "R162W"): ("TNNI3 R162W", 10, -16),
       ("msh2", "N127S"): ("MSH2 N127S", -8, 16),
       ("mlh1", "V384D"): ("MLH1 V384D", -10, 14),
       ("mlh1", "K618E"): ("MLH1 K618E", 10, -16),
       ("hbb", "E6V"): ("HBB E6V", 10, 10),
       ("calm1", "N98S"): ("CALM1 N98S", 6, -14),
       ("cfh", "R78G"): ("CFH R78G\n(fold-neutral\ninterface)", -10, -20)}
for _, r in pv.iterrows():
    k = (r["gene"].lower(), r["variant"])
    if k in div:
        lbl, dx, dy = div[k]
        col = "#8a2846" if k == ("cfh", "R78G") else "0.15"
        ax.annotate(lbl, (r["AM_score"], r["_y"]), fontsize=6.6,
                    xytext=(dx, dy), textcoords="offset points",
                    ha="left" if dx > 0 else "right", color=col,
                    arrowprops=dict(arrowstyle="-", color="0.55", lw=0.7), zorder=5)

ax.set_yticks([1, 2, 3, 4])
ax.set_yticklabels(["T1\nstrong", "T2\nmoderate", "T3\nuncertain", "T4\nminimal"])
ax.set_ylim(4.6, 0.4)
ax.set_xlim(0, 1.02)
ax.set_xlabel("AlphaMissense pathogenicity score")
ax.set_ylabel("COMAVI structural-evidence tier")
ax.legend(loc="lower right", frameon=False, fontsize=7.5)
ax.text(0.01, 0.735, "VWF R1334Q, A1381T: no AM score\n"
        "(protein >2700 aa, outside AM release)",
        ha="left", va="top", fontsize=5.8, color="0.5", style="italic")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"plotted n={len(pv)} (pathogenic {npath}, benign {nben})")
print(f"saved {OUT.relative_to(REPO)}")
