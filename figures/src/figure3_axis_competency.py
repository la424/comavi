#!/usr/bin/env python3
"""COMAVI Figure 3 — per-axis and per-mechanism-class competency across the
FoldX calling-threshold sweep (61-variant canonical benchmark).

Panel a: structural agreement per axis (monomer-fold / complex-fold / binding).
Panel b: mechanism-consistency per mechanism class.

Panel a numerators/denominators are the frozen canonical values recorded in
docs/COMAVI_v7_canonical_benchmark_ledger.md; this script renders them and does
not recompute that sweep (the sweep itself is run by the pipeline, see run.py).
Panel b IS recomputed here from reference_outputs/scored_61var_canonical.csv and
supplement/brct/brct_foldx_concordance.csv, and asserts that the per-class graded
denominators sum to the headline graded set.

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
#
# v7.3: recomputed from the canonical table via the pipeline's own per-axis
# decomposition rather than transcribed from the ledger. The previous hardcoded
# block silently went stale when the BRCA1-BRCT cohort was pooled into the
# graded arms (monomer denominator 16 -> 27). Deriving it here means the panel
# cannot disagree with the headline it decomposes.
import sys
import pandas as pd

sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as _ac

_CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
_d = pd.read_csv(_CANON, low_memory=False)
_partners = [p for p in _ac.discover_partners(_d) if p]
_TAG2COL = dict(zip(TAGS, ["t10", "t15", "t20", "t25", "tSAP"]))


def _thr_triple(spec):
    if isinstance(spec, dict):
        return spec["monomer"], spec["fold"], spec["binding"]
    return spec, spec, spec


_SPECS = dict(_ac.THRESHOLD_SPECS)
SWEEP = {}
for _t in TAGS:
    _acc = {k: [0, 0] for k in ("monomer", "fold", "binding", "tier")}
    for _, _r in _d.iterrows():
        _per = _ac.structural_agreement_by_axis(
            _r, _partners, *_thr_triple(_SPECS[_TAG2COL[_t]]))
        for _k, (_a, _b) in _per.items():
            _acc[_k][0] += _a
            _acc[_k][1] += _b
    SWEEP[_t] = {k: tuple(v) for k, v in _acc.items()}

# The panel plots the three DDG axes; the tier axis is the fourth term of the
# headline and is shown separately (Fig. 2), so assert the four close to it.
_hd = sum(SWEEP["t2.5"][k][1] for k in ("monomer", "fold", "binding", "tier"))
_hn = sum(SWEEP["t2.5"][k][0] for k in ("monomer", "fold", "binding", "tier"))
assert (_hn, _hd) == (99, 131), f"panel a decomposition {_hn}/{_hd} != headline 99/131"

_ND = {k: SWEEP["t2.5"][k][1] for k in ("monomer", "fold", "binding")}
AXINFO = [("monomer", f"Monomer-fold (n={_ND['monomer']})", C_MONO),
          ("fold", f"Complex-fold (n={_ND['fold']})", C_CX),
          ("binding", f"Binding / PPI (n={_ND['binding']})", C_BIND)]

# Panel b — mechanism-consistency by mechanism class, same threshold sweep.
# Computed from the canonical table so the figure cannot drift from the manuscript.
# Denominators are GRADED variants only: abstentions are excluded, not scored as
# failures (identical convention to the headline metric), so the four class
# denominators sum exactly to the headline graded set (v7.3: n = 57). The
# BRCA1-BRCT cohort is no longer a separate reference arm — it is pooled into
# the graded headline and appears within these classes (ledger sec.17).
_BRCT = REPO / "supplement" / "brct" / "brct_foldx_concordance.csv"
_MCMAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
_COLTAGS = ["t10", "t15", "t20", "t25", "tSAP"]
_SILENT = ["structurally_silent", "structurally_uncommitted",
           "interface_uncommitted_magnitude"]


def _mc(mask):
    return [round(_d.loc[mask, f"mech_consistency_{t}"].map(_MCMAP).mean(), 4)
            for t in _COLTAGS]


def _ngraded(mask):
    return int(_d.loc[mask, "mech_consistency_t25"].map(_MCMAP).notna().sum())


# Panel b classes are keyed on `axis_signature` — the direction-explicit,
# monomer/complex-distinguished restatement of the same ground-truth axis
# pattern that `expected_mech_class` encodes. The two agree row-for-row (the
# canonical table carries both); `axis_signature` is used here because the
# legacy labels collapse fold direction and merge the two fold axes.
_CLASSES = [
    ("No structural effect", "#7a7a7a",
     _d["axis_signature"].isin(["no_structural_effect", "uncommitted"])),
    ("Complex-fold + binding destab.", "#8e6fb0",
     _d["axis_signature"].eq("complex_fold_and_binding_destab")),
    ("Binding destab., fold intact", "#c26a2a",
     _d["axis_signature"].eq("binding_destab_fold_intact")),
    ("Fold destab. (monomer+complex)", "#4a6b8a",
     _d["axis_signature"].eq("fold_destab_monomer_and_complex")
     & _d["expected_mech_class"].notna()),
]
CLASSB = {f"{lbl} (n={_ngraded(m)})": (col, _mc(m)) for lbl, col, m in _CLASSES}

# Panel c (monomer-fold direction agreement vs measured dG_U-F) was CUT at v21.
# It binarized the same 10 measured destabilizers that Figure 4a plots as a proper
# correlation against measured dG_U-F (rho=0.72), and its 4 fold-intact controls
# are already drawn in Figure 4a under their own marker styles (M1663K, P1806A as
# measured neutral controls; R1699L/Q as fold-intact/function-lost). It therefore
# added no data, while its agreement line sat below panel a's monomer-fold line
# for purely definitional reasons (binary one-axis direction vs direction-aware
# per-axis agreement), manufacturing a reconciliation trap for readers.

_tot = sum(_ngraded(m) for _, _, m in _CLASSES)
_head = int(_d["mech_consistency_t25"].map(_MCMAP).notna().sum())
assert _tot == _head, f"class denominators {_tot} != headline graded {_head}"

apply_figure_style()
fig = plt.figure(figsize=(9.2, 4.5))
gs = fig.add_gridspec(1, 2, wspace=0.26, left=0.078, right=0.985, top=0.80, bottom=0.17)

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
    if vals[4] is not None:
        ax2.plot([XPOS[3], XPOS[4]], [vals[3], vals[4]], "--", color=col, lw=1.2,
                 alpha=0.8, zorder=2)
        ax2.plot(XPOS[4], vals[4], marker="*", color=col, ms=11, zorder=4,
                 mec="white", mew=0.4)
ax2.axvline(2.5, color="#bbb", lw=1, ls=":", zorder=1)
ax2.set_xticks(XPOS); ax2.set_xticklabels(XLAB)
ax2.set_xlim(0.8, 3.55); ax2.set_ylim(0.15, 1.06)
ax2.set_xlabel("FoldX |\u0394\u0394G| call threshold (kcal/mol)")
ax2.set_ylabel("Mechanism-consistency")
ax2.set_title("Multi-axis mechanism-consistency by axis signature", loc="left")
ax2.text(2.55, 0.185, "canonical t=2.5", fontsize=6.5, color="#999", ha="left", va="bottom")
h2, l2 = ax2.get_legend_handles_labels()
h2 += [Line2D([], [], marker="*", color="#666", ms=9, lw=0)]
l2 += ["SAP = Sapozhnikov CI"]
ax2.legend(h2, l2, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=6.3,
           handlelength=1.4, labelspacing=0.3, borderpad=0.4)
set_frame(ax2); panel_letter(ax2, "b")

fig.suptitle("COMAVI benchmark competency (61-variant set; direction-aware "
             "agreement reproduces shipped-pipeline metric)",
             x=0.078, ha="left", fontsize=9.6, y=0.955)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
print("panel a t=2.5: monomer %d/%d, fold %d/%d, binding %d/%d"
      % (*SWEEP["t2.5"]["monomer"], *SWEEP["t2.5"]["fold"], *SWEEP["t2.5"]["binding"]))
print(f"panel b (computed; class graded sum {_tot} == headline {_head}):")
for _lbl, (_c, _v) in CLASSB.items():
    print(f"    {_lbl:34s} {_v}")
print(f"saved {OUT.relative_to(REPO)}")
