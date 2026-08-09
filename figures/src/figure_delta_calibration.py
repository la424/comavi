#!/usr/bin/env python3
"""COMAVI Figure — FoldX vs directly measured DDG.

Panel a: measured vs predicted, identity + fitted line, threshold crosshairs.
Panel b: signed error vs size of the true effect.
Panel c: |error| across calls-agree / borderline / substantive disagreement.

Only two of the three COMAVI axes have measured comparators. There is NO
directly measured complex-fold comparator anywhere in this benchmark.

The Glu73 barnase electrostatics cluster (g73 == True) is drawn as open red
diamonds and EXCLUDED from every fitted statistic (slope, rho, medians).

Data:  reference_outputs/COMAVI_delta_calibration_points.csv
Stats: reference_outputs/COMAVI_delta_calibration_stats.json
Run from repo root:  python figures/src/figure_delta_calibration.py
"""
import json, pathlib
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import sys as _sys
_sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _comavi_style import (COLUMN_W as CW, BASE, SECOND, TICK,
                           apply_figure_style, check_legibility, check_overlaps)

REPO = pathlib.Path(__file__).resolve().parents[2]

mpl.rcParams.update({
    'axes.labelsize': 9.0,
    'axes.linewidth': 0.6,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.titlelocation': 'left',
    'axes.titlesize': 9.0,
    'font.size': BASE,
    'legend.fontsize': 7.5,
    'legend.frameon': False,
    'lines.linewidth': 1.2,
    'patch.linewidth': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'savefig.bbox': 'tight',
    'savefig.dpi': 300.0,
    'xtick.labelsize': 7.0,
    'xtick.major.size': 3.0,
    'xtick.major.width': 0.6,
    'ytick.labelsize': 7.0,
    'ytick.major.size': 3.0,
    'ytick.major.width': 0.6
})


def set_frame(ax, style="open"):
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, fontsize=None):
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    ax.text(dx, dy, letter.lower(), transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


D = pd.read_csv(REPO / "reference_outputs/COMAVI_delta_calibration_points.csv")
STATS = json.load(open(REPO / "reference_outputs/COMAVI_delta_calibration_stats.json"))
FIT = D[~D.g73]
sl, ic = STATS["slope"], STATS["intercept"]
rho, prho = STATS["rho"], STATS["p"]
rho_delta = STATS["rho_delta_vs_measured"]
GREY = "#8a8a8a"; C_FOLD = "#d95f02"; C_BIND = "#1f6fb4"; C_EXC = "#c0392b"
MK = {"BRCA1 BRCT": "o", "Hb tetramer": "o", "Barnase–barstar": "^", "TEM1–BLIP": "s"}
COL = {"Monomer fold": C_FOLD, "Binding": C_BIND}

def _pa(ax):
    lo,hi=-1.6,9.6
    ax.add_patch(Rectangle((2.5,lo),hi-2.5,2.5-lo,facecolor="#f2f2f2",edgecolor="none",zorder=0))
    ax.add_patch(Rectangle((lo,2.5),2.5-lo,hi-2.5,facecolor="#f2f2f2",edgecolor="none",zorder=0))
    ax.plot([lo,hi],[lo,hi],ls=":",lw=0.9,color=GREY,zorder=1)
    xs=np.array([lo,hi]); ax.plot(xs,sl*xs+ic,lw=1.5,color="#333",zorder=3)
    ax.axhline(2.5,lw=0.7,color=GREY,ls="--",zorder=1); ax.axvline(2.5,lw=0.7,color=GREY,ls="--",zorder=1)
    for (ax_,sy),g in D[~D.g73].groupby(["axis","system"]):
        ax.scatter(g.measured_kcal,g.foldx_ddg,s=30,marker=MK[sy],facecolor=COL[ax_],
                   edgecolor="white",linewidth=0.5,alpha=0.92,zorder=4)
    gx=D[D.g73]
    ax.scatter(gx.measured_kcal,gx.foldx_ddg,s=34,marker="D",facecolor="none",
               edgecolor=C_EXC,linewidth=1.1,zorder=5)
    ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
    # shared identity range means the corner x/y tick labels coincide;
    # suppress the lowest x tick so they do not overprint
    ax.set_xticks([t for t in ax.get_xticks() if t > lo + 0.35*(hi-lo)/10])
    ax.set_xlabel("Measured ΔΔG (kcal/mol)"); ax.set_ylabel("FoldX ΔΔG (kcal/mol)")
    ax.set_title("FoldX ranks with measurement\nbut compresses its scale",loc="left")
    ax.text(hi-0.25,lo+0.85,f"slope {sl:.2f}   ρ = {rho:.2f}",ha="right",va="bottom",
            fontsize=SECOND,color="#333",
            bbox=dict(boxstyle="square,pad=0.15",facecolor="white",edgecolor="none",alpha=0.85))
    ax.annotate("identity",xy=(8.2,8.2),xytext=(6.0,9.3),fontsize=SECOND,color=GREY,ha="center",
                arrowprops=dict(arrowstyle="-",lw=0.6,color=GREY))
    h=[Line2D([],[],ls="none",marker=MK[s],color=COL["Monomer fold" if s=="BRCA1 BRCT" else "Binding"],
       markersize=5.5,markeredgecolor="white",markeredgewidth=0.5,
       label=s+(" (fold)" if s=="BRCA1 BRCT" else " (binding)"))
       for s in ["BRCA1 BRCT","Hb tetramer","Barnase–barstar","TEM1–BLIP"]]
    h+=[Line2D([],[],ls="none",marker="D",markerfacecolor="none",markeredgecolor=C_EXC,
               markersize=5.5,label="Glu73 cluster (excluded)"),
        Patch(facecolor="#f2f2f2",edgecolor="none",label="calls disagree at 2.5")]
    ax.legend(handles=h,loc="upper left",bbox_to_anchor=(-0.012,1.015),frameon=False,fontsize=SECOND,
              handletextpad=0.35,borderpad=0.15,labelspacing=0.3)
    set_frame(ax); panel_letter(ax,"a")


def panel_a(ax):
    _pa(ax)
    lg=ax.get_legend(); lg.set_bbox_to_anchor((0.075,1.015),transform=ax.transAxes)


def panel_b(ax):
    sub=D[~D.g73]
    for (ax_,sy),g in sub.groupby(["axis","system"]):
        ax.scatter(g.measured_kcal,g.delta,s=30,marker=MK[sy],facecolor=COL[ax_],
                   edgecolor="white",linewidth=0.5,alpha=0.92,zorder=4)
    gx=D[D.g73]
    ax.scatter(gx.measured_kcal,gx.delta,s=34,marker="D",facecolor="none",
               edgecolor=C_EXC,linewidth=1.1,zorder=5)
    ax.axhline(0,lw=0.9,color="#333",zorder=2)
    xs=np.linspace(-1.6,9.6,50); ax.plot(xs,(sl-1)*xs+ic,lw=1.2,color=GREY,zorder=3)
    ax.set_xlim(-1.6,9.6); ax.set_ylim(-6.6,3.4)
    ax.set_xlabel("Measured ΔΔG (kcal/mol)"); ax.set_ylabel("FoldX − measured (kcal/mol)")
    ax.set_title("Error grows with the size of\nthe true effect",loc="left")
    ax.text(9.4,-5.55,f"ρ = {rho_delta:.2f}",ha="right",va="bottom",fontsize=SECOND,color="#333")
    ax.text(-1.3,3.15,"FoldX over-predicts",fontsize=SECOND,color="#666",va="top")
    ax.text(-1.3,-6.35,"FoldX under-predicts",fontsize=SECOND,color="#666",va="bottom",ha="left")
    ax.annotate("expected from\nslope 0.43",xy=(7.2,(sl-1)*7.2+ic),xytext=(6.6,-1.15),fontsize=SECOND,
                color=GREY,ha="center",arrowprops=dict(arrowstyle="-",lw=0.6,color=GREY))
    set_frame(ax); panel_letter(ax,"b")


def panel_c(ax):
    sub=FIT.copy()
    cats=[("Agree",sub[~sub.straddles],"#b9b9b9"),
          ("Borderline",sub[sub.straddles&sub.both_near],"#2c7fb8"),
          ("Substantive",sub[sub.straddles&~sub.both_near],"#d7301f")]
    for i,(nm,g,c) in enumerate(cats):
        x=np.random.default_rng(0).normal(i,0.075,len(g))
        ax.scatter(x,g.absdelta,s=26,facecolor=c,edgecolor="white",linewidth=0.5,alpha=0.92,zorder=3)
        m=g.absdelta.median()
        ax.plot([i-0.26,i+0.26],[m,m],lw=2.0,color="#222",zorder=4,solid_capstyle="butt")
        ax.text(i,-0.45,f"n={len(g)}",ha="center",fontsize=SECOND,color="#444")
        ax.text(i+0.31,m,f"{m:.2f}",ha="left",va="center",fontsize=SECOND,color="#222",zorder=5)
    ax.axhline(1.0,ls="--",lw=0.7,color=GREY,zorder=1)
    ax.text(-0.55,1.08,"1 kcal/mol",ha="left",va="bottom",fontsize=SECOND,color="#666")
    ax.set_xticks(range(3)); ax.set_xticklabels([c[0] for c in cats],fontsize=SECOND)
    ax.set_xlim(-0.62,2.95); ax.set_ylim(-0.75,7.0)
    ax.set_ylabel("|FoldX − measured| (kcal/mol)")
    ax.set_title("Only borderline disagreements are\nnear-misses in energy",loc="left")
    set_frame(ax); panel_letter(ax,"c",dx=-0.24)


fig=plt.figure(figsize=(CW,3.05))
gs=fig.add_gridspec(1,3,width_ratios=[1.06,0.92,1.24],wspace=0.46,
                    left=0.072,right=0.995,top=0.80,bottom=0.185)
axA,axB,axC=[fig.add_subplot(gs[0,i]) for i in range(3)]
panel_a(axA); panel_b(axB); panel_c(axC)
fig.savefig("figures/COMAVI_Figure_delta_calibration.png",dpi=300)
fig.savefig("figures/COMAVI_Figure_delta_calibration.pdf")
print("below-floor text:", check_legibility(fig))
print("overlaps:", check_overlaps(fig))
