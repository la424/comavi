#!/usr/bin/env python3
"""Submission figures for the COMAVI manuscript, drawn to PLOS specification.

WHY THIS EXISTS
---------------
Until v7.14 the manuscript's eleven figures came from five generators in five
styles. Most were drawn wider than the text column, so the word processor
shrank them and nearly every label reached the page below the journal's 8 pt
minimum, in DejaVu Sans or Helvetica instead of an allowed font. This script
draws the five main and five supporting figures once, at the width they are
printed, in Arial, with every text element between 8 and 12 pt at that size.
It refuses to save a figure that breaks those rules, that carries a title
inside the image, or whose labels collide.

Every plotted value is read from a committed record. Nothing is typed here
except the AlphaMissense ambiguous-range limits (0.34 and 0.564), which are
constants of that method, and the Fig 1 layout, which is the author's slide
geometry stored in fig1_layout.json.

OUTPUTS (figures/submission/)
    Fig1.tif ... Fig5.tif          600 dpi, RGB, LZW: the files uploaded to PLOS
    Fig1.png ... Fig5.png          300 dpi: the copies embedded in the manuscript
    S1_Fig.tif ... S5_Fig.tif      supporting figures, same standard (+ .png)
    SHA256SUMS.txt

Run from the repo root:
    python figures/submission/make_submission_figures.py            # write everything
    python figures/submission/make_submission_figures.py --check    # rebuild in memory, compare bytes
    python figures/submission/make_submission_figures.py --only Fig3

Byte comparison needs the same Arial, FreeType and matplotlib as the machine
that wrote the files, so CI checks the committed files with
scripts/verify_figure_specs.py instead of regenerating them.
"""
import argparse
import hashlib
import io
import json
import pathlib
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text as mtext  # noqa: E402
import matplotlib.transforms as mtrans  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
from PIL import Image  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[2]
RO = REPO / "reference_outputs"
AN = RO / "isds_v1"
OUT = REPO / "figures" / "submission"
LAYOUT = OUT / "fig1_layout.json"

TIFF_DPI, PNG_DPI = 600, 300
PT_MIN, PT_MAX = 8.0, 12.0
W_MIN, W_MAX, H_MAX = 2.63, 7.5, 8.75
FULL, HALF = 6.5, 3.4
PAD = 0.03                       # inches of white kept around the cropped figure
AM_AMBIGUOUS = (0.34, 0.564)     # AlphaMissense ambiguous range (Cheng et al. 2023)

INK, GREY, RULE, BAND = "#222222", "#6E6E6E", "#B3B3B3", "#E3E3E3"
STRUCT, NOLES = "#B2182B", "#5F5F5F"              # structural mechanism / no lesion
MONO, FOLD, BIND, TIER = "#3C78D8", "#6AA84F", "#E69138", "#674EA7"   # Fig 1 colours
ENERGY = "#0E7C7B"
TAGS = ["t10", "t15", "t20", "t25", "tSAP"]
# The fifth threshold is the per-axis 95% prediction-interval upper limit. It
# was labelled "Bound", which is ambiguous in this paper specifically: the
# Introduction says "Most proteins work while bound to partners", so a reader
# meets "Bound" on an energy axis having just been taught the bound state. The
# label is now the per-axis values themselves, which also makes the axis
# homogeneous -- every tick is a threshold in kcal/mol -- and matches the row
# label the S3 Table threshold_sweep sheet already uses.
XLAB = ["1.0", "1.5", "2.0", "2.5", "2.9/2.9/3.5"]
GRADE = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
SYSNAME = {"brca1_brct": "BRCA1 BRCT", "brca1_bard1": "BRCA1–BARD1", "cam_cav12": "Calmodulin–CaV1.2",
           "cfh_c3b": "CFH–C3b", "hemoglobin_dimer": "Hemoglobin α–β dimer",
           "hemoglobin_tetramer": "Hemoglobin tetramer", "kras_craf": "KRAS–RAF1", "mlh1_pms2": "MLH1–PMS2",
           "msh2_msh6": "MSH2–MSH6", "pi3k": "PI3K p110α–p85α", "smad4_smad3": "SMAD4–SMAD3",
           "troponin_ic": "Troponin I–troponin C", "vhl_elonginc": "VHL–elongin C", "vwf_gpiba": "VWF A1–GPIbα"}


def style():
    plt.rcdefaults()
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "legend.fontsize": 8, "axes.titlesize": 8, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
        "ytick.major.size": 2.5, "xtick.major.pad": 2, "ytick.major.pad": 2, "lines.linewidth": 1.2,
        "lines.markersize": 4, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": True, "legend.facecolor": "white", "legend.edgecolor": "none",
        "legend.framealpha": 1.0, "legend.fancybox": False, "legend.borderpad": 0.2,
        "legend.handlelength": 1.8, "legend.borderaxespad": 0.3,
        "axes.unicode_minus": True, "axes.formatter.use_mathtext": False, "figure.dpi": 100,
        "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "text.color": INK, "axes.labelpad": 3, "hatch.linewidth": 0.5})


def J(path):
    return json.loads(pathlib.Path(path).read_text())


def letter(ax, s, dx=-0.36, dy=0.05):
    """Panel letter, offset in inches from the axes' top-left corner."""
    tr = ax.transAxes + mtrans.ScaledTranslation(dx, dy, ax.figure.dpi_scale_trans)
    ax.text(0, 1, s, transform=tr, fontsize=11, fontweight="bold", va="bottom", ha="left")


def ref_line(ax):
    ax.axvline(3, color=RULE, lw=0.7, ls=(0, (3, 2)), zorder=0)


def thresh_axis(ax, xlabel=True):
    ax.set_xticks(range(5), XLAB)
    ax.set_xlim(-0.35, 4.35)
    if xlabel:
        ax.set_xlabel("Energy threshold (kcal/mol)")
    ref_line(ax)


def beeswarm(ax, values, size_pt2, gap=1.08, vertical=False):
    """Offsets (data units, across the value axis) so that no two markers overlap on the page.

    The value axis is x unless vertical=True. Axis limits and layout must be final."""
    fig = ax.figure
    bb = ax.get_position()
    w_in, h_in = bb.width * fig.get_figwidth(), bb.height * fig.get_figheight()
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    sx, sy = w_in / abs(x1 - x0), h_in / abs(y1 - y0)
    if vertical:
        sx, sy = sy, sx
    diam = gap * np.sqrt(size_pt2) / 72.0
    placed, off = [], np.zeros(len(values))
    for i, v in enumerate(values):
        for k in range(200):
            cand = ((k + 1) // 2) * (1 if k % 2 else -1) * diam * 0.55
            if all((sx * (v - pv_)) ** 2 + (cand - po) ** 2 >= diam ** 2 for pv_, po in placed):
                break
        placed.append((v, cand))
        off[i] = cand / sy
    return off


def roc(y, s):
    order = np.argsort(-s, kind="mergesort")
    s_, y_ = s[order], y[order]
    idx = np.r_[np.where(np.diff(s_))[0], len(y_) - 1]
    tps = np.cumsum(y_)[idx]
    fps = (idx + 1) - tps
    return np.r_[0, fps / (len(y) - y.sum())], np.r_[0, tps / y.sum()]


def canonical_grades():
    can = pd.read_csv(RO / "scored_61var_canonical.csv")
    g = can[can.mech_consistency_t25.notna()]
    return g, g.expected_mech_class != "structurally_silent"


# ------------------------------------------------------------------ Fig 1

# The slide's box labels, mapped to the terms the manuscript uses. The geometry
# in fig1_layout.json stays a faithful extraction of the author's slide; only
# the words change, and a label missing from this table stops the build, so the
# figure cannot drift from the text if the layout is re-extracted.
FIG1_TERMS = {
    "Missense Variant": "Missense variant",
    "Isolated/Assembled Structural Models": "Isolated and assembled structural models",
    "Stability in Monomer": "Isolated-subunit stability",
    "Stability in Multimer": "Assembled-complex stability",
    "Binding Energy by Partner": "Binding energy for each partner",
    "Structural Context Tier": "Structural-context tier",
    "Priority Score (ISDS)": "Priority score",
    "Mechanism Profile": "Mechanism profile",
    "Prioritization": "Prioritization",
    "Experimental Design": "Experimental design",
    "Input": "Input",
    "Structural Evidence": "Structural evidence",
    "COMAVI Outputs": "COMAVI outputs",
    "Practical Use": "Practical use",
    "Ranks strength of modeled structural-disruption evidence": "Ranks strength of modeled structural-disruption evidence",
    "Assigns per-variant mechanisms and context": "Assigns per-variant mechanisms and context",
    "Choose variants for further evaluation": "Choose variants for further evaluation",
    "Test stability, assembly, and/or interaction": "Test stability, assembly or binding",
}


def term(text):
    if text not in FIG1_TERMS:
        raise SystemExit("Fig 1 label %r has no entry in FIG1_TERMS" % text)
    return FIG1_TERMS[text]


def fig1():
    L = J(LAYOUT)
    boxes, labels, conns = L["boxes"], L["labels"], L["connectors"]
    hdr = [lab for lab in labels if lab["box"] is None]
    x0 = min(b["x"] for b in boxes) - 0.03
    x1 = max(b["x"] + b["w"] for b in boxes) + 0.03
    y0 = min(lab["y"] for lab in hdr) + lab_inset(hdr[0]) - 0.04
    y1 = max(b["y"] + b["h"] for b in boxes) + 0.03
    S = FULL / (x1 - x0)
    fig = plt.figure(figsize=(FULL, (y1 - y0) * S))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.axis("off")

    def pt(slide_pt):
        # Scale with the drawing, round down to a half point, never below 8 pt.
        return max(PT_MIN, np.floor(slide_pt * S * 2) / 2)

    for b in boxes:
        ax.add_patch(FancyBboxPatch(
            (b["x"], b["y"]), b["w"], b["h"],
            boxstyle="round,pad=0,rounding_size=%.4f" % (b["corner"] * min(b["w"], b["h"])),
            fc=b["fill"], ec=b["line"], lw=b["lw_pt"] * S, zorder=1))
    for c in conns:
        (xa, ya), (xb, yb) = c["start"], c["end"]
        if c["arrow_at_end"]:
            ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb), arrowstyle="-|>", mutation_scale=5.5,
                                         lw=c["lw_pt"] * S, color=c["color"], shrinkA=0, shrinkB=0,
                                         zorder=2))
        else:
            ax.add_line(Line2D([xa, xb], [ya, yb], lw=c["lw_pt"] * S, color=c["color"],
                               solid_capstyle="butt", zorder=2))
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    for h in hdr:
        cx = h["x"] + h["w"] / 2
        cy = h["y"] + lab_inset(h) + 1.2 * h["pt"] / 72 / 2
        ax.text(cx, cy, term(h["text"]), ha="center", va="center", fontsize=pt(h["pt"]), fontweight="bold",
                zorder=3)
    for b in boxes:
        own = sorted([lab for lab in labels if lab["box"] == b["id"]], key=lambda d: d["y"])
        title, subs = own[0], own[1:]
        blocks = []
        for lab in [title] + subs:
            size = pt(lab["pt"])
            avail = min(lab["w"] - 2 * lab_inset(lab), b["w"] - 2 * 0.05) * S
            wrapped = wrap(fig, ren, term(lab["text"]), size, avail)
            nlines = wrapped.count("\n") + 1
            blocks.append((wrapped, size, nlines * size * 1.18 / 72))
        gap = 0.035
        total = sum(h for _, _, h in blocks) + gap * (len(blocks) - 1)
        y = b["y"] + b["h"] / 2 - total / S / 2
        for wrapped, size, h in blocks:
            ax.text(b["x"] + b["w"] / 2, y, wrapped, ha="center", va="top", fontsize=size,
                    fontweight="bold", linespacing=1.18, zorder=3)
            y += (h + gap) / S
    return fig


def lab_inset(lab):
    return lab.get("inset_in", 0.1)


def wrap(fig, ren, text, size, max_w_in):
    prop = font_manager.FontProperties(family="Arial", weight="bold", size=size)
    lines, cur = [], ""
    # Break at spaces, and after a slash when a word alone is wider than the box.
    words = []
    for word in text.split():
        w_px, _, _ = ren.get_text_width_height_descent(word, prop, ismath=False)
        if "/" in word and w_px / fig.dpi > max_w_in:
            head, tail = word.split("/", 1)
            words += [head + "/", tail]
        else:
            words.append(word)
    for word in words:
        trial = (cur + ("" if cur.endswith("/") else " ") + word).strip()
        w, _, _ = ren.get_text_width_height_descent(trial, prop, ismath=False)
        if w / fig.dpi <= max_w_in or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = word
    lines.append(cur)
    return "\n".join(lines)


# ------------------------------------------------------------------ priority score (Fig 4)

def fig_priority():
    pv = pd.read_csv(AN / "ISDS_v1_per_variant.csv")
    tk = pd.read_csv(AN / "ISDS_v1_top_k.csv")
    pop = J(AN / "ISDS_v1_summary.json")["population"]
    y = pv.structural_ground_truth.astype(bool)
    assert (int(y.sum()), int((~y).sum())) == (pop["positive"], pop["negative"])
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3), gridspec_kw=dict(width_ratios=[1.12, 1]))
    fig.subplots_adjust(left=0.165, right=0.99, bottom=0.2, top=0.93, wspace=0.42)
    a.set_xlim(0, 1)
    a.set_ylim(1.5, -0.5)
    for row, mask, col in [(0, y, STRUCT), (1, ~y, NOLES)]:
        v = np.sort(pv.loc[mask, "isds_v1"].to_numpy())
        a.scatter(v, row + beeswarm(a, v, 13), s=13, color=col, edgecolor="white", linewidth=0.4,
                  zorder=3)
    a.set_yticks([0, 1], ["Structural\nmechanism (n = %d)" % pop["positive"],
                          "No lesion\n(n = %d)" % pop["negative"]])
    a.set_xlabel("Priority score")
    a.tick_params(axis="y", length=0)
    a.spines["left"].set_visible(False)
    letter(a, "A", dx=-0.95)
    for score, col, ls, mk, lab in [("isds_v1", INK, "-", "o", "Priority score"),
                                    ("isds_context_component", TIER, "--", "s", "Context component"),
                                    ("isds_energy_component", ENERGY, ":", "^", "Energy component")]:
        t = tk[tk.score == score].sort_values("k")
        assert (t.k == t.k_effective).all()
        b.plot(t.k, t.structural_mechanisms_in_top_k, ls=ls, color=col, marker=mk, ms=3.6,
               mfc=col, mec=col, lw=1.2, label=lab, zorder=3)
    ks = sorted(tk.k.unique())
    b.set_xticks(ks)
    b.set_xlim(ks[0] - 1.5, ks[-1] + 1.5)
    b.set_ylim(0, pop["positive"] + 0.9)
    b.set_yticks([0, 5, 10, 15, pop["positive"]])
    b.set_xlabel("Rank depth (top k variants)")
    b.set_ylabel("Structural mechanisms\nrecovered")
    b.legend(loc="lower right")
    letter(b, "B", dx=-0.52)
    fig.subplots_adjust(left=0.165, right=0.99, bottom=0.2, top=0.93, wspace=0.42)
    return fig


# ------------------------------------------------------------------ threshold sweep (Fig 2)

def fig_threshold():
    sw = J(RO / "COMAVI_threshold_sweep.json")
    rows = sw["rows"]
    assert [r["threshold_tag"] for r in rows] == TAGS
    g, struct = canonical_grades()
    arm_s = [g.loc[struct, "mech_consistency_%s" % t].map(GRADE).mean() for t in TAGS]
    arm_n = [g.loc[~struct, "mech_consistency_%s" % t].map(GRADE).mean() for t in TAGS]
    op = pd.read_csv(RO / "COMAVI_threshold_operating_points.csv")
    assert op.threshold.tolist() == TAGS
    assert np.allclose(arm_s, op.sensitivity, atol=6e-4) and np.allclose(arm_n, op.specificity, atol=6e-4)
    me = J(RO / "COMAVI_measured_effect_tables.json")["measured_effects_recovered"]["by_threshold"]
    assert [m["threshold"] for m in me] == ["1.0", "1.5", "2.0", "2.5", "tSAP"]
    X = np.arange(5)
    fig, ax = plt.subplots(2, 2, figsize=(FULL, 4.25))
    # Panels are laid out in CITATION order. The Results cite the
    # mechanism-pattern score first, then the three binary rates, then per-axis
    # direction agreement, then measured effects -- but the panels used to be
    # drawn rates/agreement/score/effects, so the text referred to 2C before 2A.
    # The variable names below still carry the same content as before; only
    # which grid cell each one occupies, and the letter it is given, changed.
    # Unpacking in this order moves the score panel to top-left without
    # touching any of the drawing code.
    c, a, b, d = ax.ravel()

    for key, col, ls, mk, mfc, lab in [("detection", STRUCT, "-", "o", STRUCT, "Detection"),
                                       ("attribution", STRUCT, "--", "o", "white", "Attribution"),
                                       ("correct_rejection", NOLES, "-", "^", NOLES, "Correct rejection")]:
        a.plot(X, [r[key] for r in rows], ls=ls, color=col, marker=mk, mfc=mfc, mec=col, ms=4,
               label=lab, zorder=3)
    a.set_ylim(0, 1.02)
    a.set_ylabel("Fraction of variants")
    a.legend(loc="lower right")
    thresh_axis(a, xlabel=False)
    letter(a, "B")

    for k, col, mk, lab in [("axis_monomer", MONO, "o", "Isolated subunit"),
                            ("axis_fold", FOLD, "s", "Assembled complex"),
                            ("axis_binding", BIND, "^", "Binding")]:
        vals = np.array([r[k + "_k"] / r[k + "_n"] for r in rows])
        b.plot(X, vals, "-", color=col, marker=mk, mfc=col, mec=col, ms=3.6, label=lab, zorder=3)
        best = np.where(np.isclose(vals, vals.max()))[0]
        b.plot(best, vals[best], ls="none", marker="o", ms=8, mfc="none", mec=col, mew=0.9, zorder=4)
    b.set_ylim(0.4, 1.0)
    b.set_ylabel("Direction agreement")
    b.legend(loc="lower right", ncol=1)
    thresh_axis(b)
    letter(b, "C")

    bands = sw["cluster_bands"]
    lo = [bands[t]["lo"] for t in TAGS]
    hi = [bands[t]["hi"] for t in TAGS]
    c.fill_between(X, lo, hi, color=BAND, lw=0, zorder=1)
    c.plot(X, [r["whole_variant"] for r in rows], "-", color=INK, marker="o", ms=4, label="Pooled",
           zorder=4)
    c.plot(X, arm_s, "--", color=STRUCT, marker="s", ms=3.6, label="Structural arm", zorder=3)
    c.plot(X, arm_n, "--", color=NOLES, marker="^", ms=3.6, label="No-lesion arm", zorder=3)
    c.set_ylim(0, 1.02)
    c.set_ylabel("Mechanism-pattern score")
    c.legend(loc="lower right")
    thresh_axis(c, xlabel=False)
    letter(c, "A")

    d.plot(X, [m["fraction"] for m in me], "-", color=INK, marker="o", ms=4, zorder=3)
    d.set_ylim(0, 1.02)
    d.set_ylabel("Measured destabilizations\nrecovered")
    thresh_axis(d)
    letter(d, "D")
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.11, top=0.95, wspace=0.34, hspace=0.36)
    return fig


# ------------------------------------------------------------------ calibration (Fig 3)

SYSTEMS = [("Barnase–barstar", "o", BIND, "Barnase–barstar"), ("TEM1–BLIP", "s", BIND, "TEM1–BLIP"),
           ("BRCA1 BRCT", "^", MONO, "BRCA1 BRCT"), ("Hb tetramer", "D", BIND, "Hemoglobin")]


def calibration_points():
    d = pd.read_csv(RO / "COMAVI_delta_calibration_points.csv")
    d["g73"] = d.g73.fillna(False).astype(bool)
    assert set(d.system) == {s[0] for s in SYSTEMS}
    return d


def draw_points(ax, d, ycol):
    for sysname, mk, col, _ in SYSTEMS:
        s = d[(d.system == sysname) & ~d.g73]
        ax.scatter(s.measured_kcal, s[ycol], marker=mk, s=15, color=col, edgecolor="white", linewidth=0.35,
                   zorder=3)
    g = d[d.g73]
    ax.scatter(g.measured_kcal, g[ycol], marker="o", s=15, facecolor="white", edgecolor=BIND,
               linewidth=0.9, zorder=3)


def fig_calibration():
    d = calibration_points()
    st = J(RO / "COMAVI_delta_calibration_stats.json")
    assert st["n_fit"] == len(d)
    fig, ax = plt.subplots(2, 2, figsize=(FULL, 4.75))
    a, b, c, e = ax.ravel()
    lo, hi = min(d.measured_kcal.min(), d.foldx_ddg.min()) - 0.5, max(d.measured_kcal.max(), d.foldx_ddg.max()) + 0.5
    for axx in (a,):
        axx.plot([lo, hi], [lo, hi], ls=(0, (4, 2)), color=RULE, lw=0.8, zorder=1)
        xs = np.array([d.measured_kcal.min(), d.measured_kcal.max()])
        axx.plot(xs, st["intercept"] + st["slope"] * xs, "-", color=INK, lw=1.1, zorder=2)
        axx.axhline(2.5, color=RULE, lw=0.6, ls=":", zorder=0)
        axx.axvline(2.5, color=RULE, lw=0.6, ls=":", zorder=0)
    draw_points(a, d, "foldx_ddg")
    a.set_xlim(lo, hi)
    a.set_ylim(lo, max(d.foldx_ddg.max(), 7) + 0.5)
    a.set_xlabel("Measured ΔΔG (kcal/mol)")
    a.set_ylabel("Predicted ΔΔG (kcal/mol)")
    letter(a, "A")

    d["err"] = d.foldx_ddg - d.measured_kcal
    b.axhline(0, color=RULE, lw=0.8, ls=(0, (4, 2)), zorder=1)
    draw_points(b, d, "err")
    b.set_xlim(lo, hi)
    b.set_xlabel("Measured ΔΔG (kcal/mol)")
    b.set_ylabel("Predicted − measured\n(kcal/mol)")
    letter(b, "B")

    br = d[d.system == "BRCA1 BRCT"]
    hb = d[d.system == "Hb tetramer"]
    for axx, sub, mk, col, xl, yl in [(c, br, "^", MONO, "Measured unfolding ΔΔG (kcal/mol)",
                                       "Predicted isolated-subunit\nΔΔG (kcal/mol)"),
                                      (e, hb, "D", BIND, "Measured assembly ΔΔG (kcal/mol)",
                                       "Predicted binding\nΔΔG (kcal/mol)")]:
        m = max(sub.measured_kcal.max(), sub.foldx_ddg.max()) + 0.6
        n = min(sub.measured_kcal.min(), sub.foldx_ddg.min(), 0) - 0.6
        axx.plot([n, m], [n, m], ls=(0, (4, 2)), color=RULE, lw=0.8, zorder=1)
        axx.axhline(2.5, color=RULE, lw=0.6, ls=":", zorder=0)
        axx.axvline(2.5, color=RULE, lw=0.6, ls=":", zorder=0)
        axx.scatter(sub.measured_kcal, sub.foldx_ddg, marker=mk, s=18, color=col, edgecolor="white",
                    linewidth=0.35, zorder=3)
        axx.set_xlim(n, m)
        axx.set_ylim(n, m)
        axx.set_xlabel(xl)
        axx.set_ylabel(yl)
    miss = br[(br.measured_kcal >= 2.5) & (br.foldx_ddg < 2.5)]
    for _, r in miss.iterrows():
        c.annotate(r.variant, (r.measured_kcal, r.foldx_ddg), xytext=(6, -2), textcoords="offset points",
                   fontsize=8, va="center")
    offs = {"N102T": (6, 0), "W37Y": (-6, 3), "W37A": (-6, 4), "W37G": (-6, 6), "W37E": (0, -10)}
    for _, r in hb.iterrows():
        dx, dy = offs[r.variant]
        e.annotate(r.variant, (r.measured_kcal, r.foldx_ddg), xytext=(dx, dy), textcoords="offset points",
                   fontsize=8, ha="left" if dx > 0 else ("right" if dx < 0 else "center"), va="center")
    letter(c, "C")
    letter(e, "D")
    handles = [Line2D([], [], ls="none", marker=mk, color=col, mec="white", mew=0.35, ms=5, label=lab)
               for _, mk, col, lab in SYSTEMS]
    handles += [Line2D([], [], ls="none", marker="o", mfc="white", mec=BIND, mew=0.9, ms=5,
                       label="Barnase Glu73"),
                Line2D([], [], color=INK, lw=1.1, label="Least-squares fit"),
                Line2D([], [], color=RULE, lw=0.8, ls=(0, (4, 2)), label="Identity")]
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.54, 1.0), columnspacing=1.2,
               handletextpad=0.4)
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.095, top=0.87, wspace=0.34, hspace=0.42)
    return fig


# ------------------------------------------------------------------ Fig 5

FIG5_LABELS = {
    "isds_v1": {"E6V": (8, -5), "R78G": (8, 0), "H1047R": (8, 0), "R162W": (6, -9), "N127S": (-8, 6)},
    "isds_energy_component": {"E6V": (8, 4), "R78G": (8, 0), "H1047R": (6, -10), "R162W": (7, -7),
                              "N127S": (-8, 6)},
}


def fig5():
    cs = pd.read_csv(AN / "ISDS_v1_alphamissense_common_set.csv")
    path = cs.clinical_y.eq(1)
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.55))
    for axx, col, ylab, lab in [(a, "isds_v1", "Priority score", "A"),
                                (b, "isds_energy_component", "Energy component", "B")]:
        axx.axvspan(*AM_AMBIGUOUS, color=BAND, lw=0, zorder=0)
        axx.scatter(cs.loc[path, "AM pathogenicity"], cs.loc[path, col], s=15, color=INK, edgecolor="white",
                    linewidth=0.35, zorder=3, label="Pathogenic")
        axx.scatter(cs.loc[~path, "AM pathogenicity"], cs.loc[~path, col], s=15, facecolor="white",
                    edgecolor=INK, linewidth=0.8, zorder=3, label="Benign")
        for v, (dx, dy) in FIG5_LABELS[col].items():
            r = cs[cs.variant == v].iloc[0]
            axx.annotate(v, (r["AM pathogenicity"], r[col]), xytext=(dx, dy), textcoords="offset points",
                         fontsize=8, ha="left" if dx > 0 else "right", va="center",
                         arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY, shrinkA=0, shrinkB=2))
        axx.set_xlim(0, 1.02)
        axx.set_ylim(0, 1.02)
        axx.set_xlabel("AlphaMissense pathogenicity")
        axx.set_ylabel(ylab)
        letter(axx, lab)
    a.legend(loc="upper left", handletextpad=0.2, borderaxespad=0.1)
    fig.subplots_adjust(left=0.085, right=0.99, bottom=0.17, top=0.94, wspace=0.3)
    return fig


# ------------------------------------------------------------------ S1 Fig

def s1():
    z = np.load(RO / "stress_tests" / "comavi_stress_draws.npz")
    st = J(RO / "stress_tests" / "COMAVI_stress_verified_statistics.json")
    loso = pd.read_csv(RO / "stress_tests" / "comavi_leave_one_system_out.csv").sort_values("mc")
    assert set(loso.dropped) == set(SYSNAME)
    obs = float(z["mc_obs"][0])
    ci = st["system_cluster_ci95"]["mc"]["full_precision"]
    assert np.allclose(np.percentile(z["cluster_mc"], [2.5, 97.5]), ci, atol=1e-9)
    W_, H_ = FULL, 4.0
    fig = plt.figure(figsize=(W_, H_))

    def box(x, y, w, h):
        return fig.add_axes([x / W_, y / H_, w / W_, h / H_])

    w = 1.325
    a, b = box(0.55, 2.3, w, 1.45), box(0.55 + w + 0.55, 2.3, w, 1.45)
    d, c = box(0.55, 0.45, 2 * w + 0.55, 1.3), box(0.55 + 2 * w + 0.55 + 1.35, 0.45, 1.3, 3.3)
    a.hist(z["mc_null"], bins=30, color="#BDBDBD", lw=0)
    a.axvline(obs, color=INK, lw=1.2)
    a.set_xlabel("Score, labels permuted")
    a.set_ylabel("Permutations")
    letter(a, "A", dx=-0.42)
    b.hist(z["cluster_mc"], bins=30, color="#BDBDBD", lw=0)
    for v in ci:
        b.axvline(v, color=INK, lw=0.8, ls=(0, (3, 2)))
    b.axvline(obs, color=INK, lw=1.2)
    b.set_xlabel("Score, systems resampled")
    b.set_ylabel("Resamples")
    letter(b, "B", dx=-0.42)
    yy = np.arange(len(loso))
    c.hlines(yy, obs, loso.mc, color="#CCCCCC", lw=0.8, zorder=1)
    c.scatter(loso.mc, yy, s=13, color=INK, zorder=3)
    c.axvline(obs, color=GREY, lw=0.8, ls=(0, (4, 2)), zorder=2)
    c.set_yticks(yy, [SYSNAME[s] for s in loso.dropped])
    c.set_ylim(-0.7, len(loso) - 0.3)
    c.set_xlabel("Score, one system left out")
    c.tick_params(axis="y", length=0)
    letter(c, "C", dx=-1.25)
    vals, cnts = np.unique(np.round(z["noise_mc"], 6), return_counts=True)
    d.vlines(vals, 0, cnts, color=GREY, lw=2.2, zorder=2)
    d.axvline(obs, color=INK, lw=1.2, zorder=1)
    d.set_xlabel("Score under replicate noise")
    d.set_ylabel("Draws")
    letter(d, "D", dx=-0.42)
    return fig


# ------------------------------------------------------------------ S2 Fig

def s2():
    bb = pd.read_csv(REPO / "supplement/skempi/bb_binding_validation.csv")
    jt = pd.read_csv(REPO / "supplement/skempi/jt_binding_validation.csv")
    rec = J(RO / "COMAVI_skempi_validation.json")
    g73 = bb.foldx_code.str.match(r"E[A-Z]?73")
    assert int(g73.sum()) == rec["glu73_cluster"]["n"]
    assert len(bb) + len(jt) == rec["combined_including_glu73"]["n"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.6))
    allv = pd.concat([bb, jt])
    lo = min(allv.ddg_meas.min(), allv.ddg_pred_bind.min()) - 0.4
    hi = max(allv.ddg_meas.max(), allv.ddg_pred_bind.max()) + 0.4
    for axx, df, mask, lab in [(a, bb, g73, "A"), (b, jt, pd.Series(False, index=jt.index), "B")]:
        axx.plot([lo, hi], [lo, hi], ls=(0, (4, 2)), color=RULE, lw=0.8, zorder=1)
        axx.axhline(0, color=RULE, lw=0.6, ls=":", zorder=0)
        axx.axvline(0, color=RULE, lw=0.6, ls=":", zorder=0)
        axx.scatter(df.loc[~mask, "ddg_meas"], df.loc[~mask, "ddg_pred_bind"], s=15, color=BIND,
                    edgecolor="white", linewidth=0.35, zorder=3)
        if mask.any():
            axx.scatter(df.loc[mask, "ddg_meas"], df.loc[mask, "ddg_pred_bind"], s=15, facecolor="white",
                        edgecolor=BIND, linewidth=0.9, zorder=3, label="Glu73")
            axx.legend(loc="upper left", handletextpad=0.2)
        axx.set_xlim(lo, hi)
        axx.set_ylim(lo, hi)
        axx.set_xlabel("Measured binding ΔΔG (kcal/mol)")
        axx.set_ylabel("Predicted binding ΔΔG (kcal/mol)")
        letter(axx, lab)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.17, top=0.94, wspace=0.32)
    return fig


# ------------------------------------------------------------------ S3 Fig

def s3():
    cb = {r["screen"]: r for r in J(AN / "ISDS_v1_summary.json")["component_binary_baselines"]}
    order = [("interface_status_alone", "Interface\nstatus alone"),
             ("tier_without_interface_bonus", "Tier without\ninterface bonus"),
             ("full_tier_1_2", "Full tier\n(Tier 1–2)")]
    fig, ax = plt.subplots(figsize=(HALF, 2.45))
    x = np.arange(len(order))
    w = 0.36
    sens = [cb[k]["sensitivity"] for k, _ in order]
    spec = [cb[k]["specificity"] for k, _ in order]
    ax.bar(x - w / 2, sens, w, color=STRUCT, label="Sensitivity")
    ax.bar(x + w / 2, spec, w, color=NOLES, label="Specificity")
    for i, (u, v) in enumerate(zip(sens, spec)):
        ax.text(i - w / 2, u + 0.02, "%.3f" % u, ha="center", va="bottom", fontsize=8)
        ax.text(i + w / 2, v + 0.02, "%.3f" % v, ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x, [lab for _, lab in order])
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Fraction")
    ax.tick_params(axis="x", length=0)
    ax.legend(loc="upper left", ncol=2, handlelength=1.0, columnspacing=1.0)
    fig.subplots_adjust(left=0.17, right=0.99, bottom=0.2, top=0.98)
    return fig


# ------------------------------------------------------------------ S4 Fig

def s4():
    pv = pd.read_csv(AN / "ISDS_v1_per_variant.csv")
    pm = pd.read_csv(AN / "ISDS_v1_primary_metrics.csv").set_index("score")
    y = pv.structural_ground_truth.astype(int).to_numpy()
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.5))
    R = np.linspace(0, 8, 400)
    a.plot(R, R / (1 + R), color=INK, lw=1.3)
    for r in (0.5, 1, 2, 4):
        a.plot([r], [r / (1 + r)], "o", color=INK, ms=3.6)
        a.annotate("R = %g" % r, (r, r / (1 + r)), xytext=(5, -9), textcoords="offset points", fontsize=8)
    a.set_xlim(0, 8)
    a.set_ylim(0, 1)
    a.set_xlabel("Largest energy ratio, R")
    a.set_ylabel("Energy component, E")
    letter(a, "A")
    b.plot([0, 1], [0, 1], ls=(0, (4, 2)), color=RULE, lw=0.8)
    for col, c_, ls, lab in [("isds_v1", INK, "-", "Priority score"),
                             ("isds_energy_component", ENERGY, ":", "Energy component"),
                             ("isds_context_component", TIER, "--", "Context component")]:
        fpr, tpr = roc(y, pv[col].to_numpy())
        auc = float(np.trapezoid(tpr, fpr)) if hasattr(np, "trapezoid") else float(np.trapz(tpr, fpr))
        assert abs(auc - pm.loc[col, "roc_auc"]) < 1e-9, (col, auc, pm.loc[col, "roc_auc"])
        b.plot(fpr, tpr, ls=ls, color=c_, lw=1.2, label="%s (AUC %.3f)" % (lab, pm.loc[col, "roc_auc"]))
    b.set_xlim(0, 1)
    b.set_ylim(0, 1.02)
    b.set_xlabel("False-positive rate")
    b.set_ylabel("True-positive rate")
    b.legend(loc="lower right", handlelength=1.6)
    letter(b, "B")
    fig.subplots_adjust(left=0.085, right=0.99, bottom=0.18, top=0.94, wspace=0.3)
    return fig


# ------------------------------------------------------------------ S5 Fig

def s5():
    d = calibration_points()
    rec = J(RO / "COMAVI_measured_effect_tables.json")["call_relationships"]
    st = d.straddles.fillna(False).astype(bool)
    bn = d.both_near.fillna(False).astype(bool)
    # The two flags are stored, not derived here. Hold them to the definitions the
    # caption states, so a changed record fails instead of regrouping silently.
    assert (st == ((d.measured_kcal >= 2.5) != (d.foldx_ddg >= 2.5))).all(), "straddles != calls differ at 2.5"
    assert (bn == (((d.measured_kcal - 2.5).abs() <= 1.0) & ((d.foldx_ddg - 2.5).abs() <= 1.0))).all(), \
        "both_near != both values within 1.0 kcal/mol of 2.5"
    err = (d.foldx_ddg - d.measured_kcal).abs()
    classes = [("Agreement", ~st), ("Borderline disagreement", st & bn), ("Substantive disagreement", st & ~bn)]
    want = {c["call_relationship"]: c for c in rec["classes"]}
    fig, ax = plt.subplots(figsize=(HALF, 2.5))
    fig.subplots_adjust(left=0.16, right=0.99, bottom=0.2, top=0.98)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(0, float(err.max()) * 1.06)
    ticks = []
    for i, (name, mask) in enumerate(classes):
        assert int(mask.sum()) == want[name]["n"]
        med = float(err[mask].median())
        assert abs(round(med, 2) - want[name]["median_abs_error"]) < 1e-9
        v = err[mask].to_numpy()
        order = np.argsort(v, kind="mergesort")
        off = np.zeros(len(v))
        off[order] = beeswarm(ax, v[order], 12, vertical=True)
        g = d.g73[mask].to_numpy()
        ax.scatter(i + off[~g], v[~g], s=12, color=INK, edgecolor="white", linewidth=0.3, zorder=3)
        ax.scatter(i + off[g], v[g], s=12, facecolor="white", edgecolor=INK, linewidth=0.8, zorder=3)
        ax.hlines(med, i - 0.32, i + 0.32, color="#9E9E9E", lw=2.6, zorder=2)
        ticks.append("%s\n(n = %d)" % (name.split()[0], int(mask.sum())))
    ax.set_xticks(range(3), ticks)
    ax.set_ylabel("|Predicted − measured| (kcal/mol)")
    ax.tick_params(axis="x", length=0)
    return fig


# Numbered in the order the Results cite them: threshold, calibration, priority; and
# for the supporting figures, call relationships, context comparators, transformation.
FIGURES = {"Fig1": fig1, "Fig2": fig_threshold, "Fig3": fig_calibration, "Fig4": fig_priority, "Fig5": fig5,
           "S1_Fig": s1, "S2_Fig": s2, "S3_Fig": s5, "S4_Fig": s3, "S5_Fig": s4}


# ------------------------------------------------------------------ specification

def font_file(t):
    return pathlib.Path(font_manager.findfont(t.get_fontproperties(), fallback_to_default=False)).name


def texts_of(fig):
    fig.canvas.draw()
    tick_ids, keep = set(), []
    for ax in fig.axes:
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            lo, hi = sorted(lim)
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                for lab in (tk.label1, tk.label2):
                    tick_ids.add(id(lab))
                    if (ax.axison and axis.get_visible() and lab.get_visible() and lab.get_text().strip()
                            and lo - 1e-9 <= tk.get_loc() <= hi + 1e-9):
                        keep.append(lab)
    other = [t for t in fig.findobj(mtext.Text)
             if id(t) not in tick_ids and t.get_visible() and t.get_text().strip()]
    return other + keep


def spec_violations(fig):
    bad = []
    if fig._suptitle is not None and fig._suptitle.get_text().strip():
        bad.append("figure title inside the image")
    for ax in fig.axes:
        for loc in ("left", "center", "right"):
            if ax.get_title(loc=loc).strip():
                bad.append("axes title inside the image: %r" % ax.get_title(loc=loc))
    ren = fig.canvas.get_renderer()
    boxes = []
    for t in texts_of(fig):
        s = t.get_text()
        if "$" in s:
            bad.append("mathtext (non-Arial glyphs): %r" % s[:40])
        size = t.get_fontsize()
        if not (PT_MIN - 1e-6 <= size <= PT_MAX + 1e-6):
            bad.append("%.2f pt outside %g-%g pt: %r" % (size, PT_MIN, PT_MAX, s[:40]))
        try:
            ff = font_file(t)
        except Exception as exc:  # noqa: BLE001
            ff = "not found (%s)" % type(exc).__name__
        if not re.match(r"(?i)^arial", ff):
            bad.append("font %s, not Arial: %r" % (ff, s[:40]))
        boxes.append((s, t.get_window_extent(ren)))
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > 1.0 and oy > 1.0:
                bad.append("labels overlap: %r / %r" % (boxes[i][0][:30], boxes[j][0][:30]))
    return bad


def render(fig, dpi):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=PAD, facecolor="white",
                metadata={"Software": None})
    im = Image.open(io.BytesIO(buf.getvalue()))
    im.load()
    return im.convert("RGB")


TIFF_TYPE_SIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8, 16: 8}


def canonical_tiff(data):
    """Zero every byte of a little-endian TIFF that no structure points to.

    libtiff pads the strip data to an even offset before the directory and
    leaves that pad byte uninitialized, so two encodes of the same pixels
    differed in exactly one byte and the files could never verify against a
    checksum. Pixels, tags and offsets are untouched; only unreferenced bytes
    are set to zero."""
    import struct
    buf = bytearray(data)
    assert buf[:4] == b"II*\x00", "expected a little-endian classic TIFF"
    used = [(0, 8)]
    ifd = struct.unpack_from("<I", buf, 4)[0]
    strips, counts = None, None
    while ifd:
        n_ent = struct.unpack_from("<H", buf, ifd)[0]
        used.append((ifd, ifd + 2 + 12 * n_ent + 4))
        for k in range(n_ent):
            tag, typ, cnt, val = struct.unpack_from("<HHII", buf, ifd + 2 + 12 * k)
            size = TIFF_TYPE_SIZE[typ] * cnt
            fmt = {3: "H", 4: "I"}.get(typ)
            if size > 4:
                used.append((val, val + size))
                vals = struct.unpack_from("<%d%s" % (cnt, fmt), buf, val) if fmt else None
            else:
                vals = struct.unpack_from("<%d%s" % (cnt, fmt), buf, ifd + 2 + 12 * k + 8) if fmt else None
            if tag == 273:
                strips = vals
            elif tag == 279:
                counts = vals
        ifd = struct.unpack_from("<I", buf, ifd + 2 + 12 * n_ent)[0]
    used += [(o_, o_ + c_) for o_, c_ in zip(strips, counts)]
    keep = bytearray(len(buf))
    for a, b in used:
        keep[a:b] = b"\x01" * (b - a)
    for i, flag in enumerate(keep):
        if not flag:
            buf[i] = 0
    return bytes(buf)


def encode(im, kind, dpi):
    out = io.BytesIO()
    if kind == "tif":
        im.save(out, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
        return canonical_tiff(out.getvalue())
    im.save(out, format="PNG", dpi=(dpi, dpi))
    return out.getvalue()


def build(name):
    style()
    fig = FIGURES[name]()
    bad = spec_violations(fig)
    if bad:
        plt.close(fig)
        raise SystemExit("%s violates the figure specification:\n  " % name + "\n  ".join(bad))
    # Render once, at the upload resolution, and derive the embed copy from it.
    # Rendering twice let the tight bounding box differ by a few pixels between
    # resolutions, so the PNG in the manuscript and the TIFF uploaded to PLOS
    # were not the same picture at the same size.
    im = render(fig, TIFF_DPI)
    plt.close(fig)
    w, h = im.size
    if w % 2 or h % 2:
        even = Image.new("RGB", (w + w % 2, h + h % 2), "white")
        even.paste(im, (0, 0))
        im = even
    w_in, h_in = im.size[0] / TIFF_DPI, im.size[1] / TIFF_DPI
    if not (W_MIN <= w_in <= W_MAX and h_in <= H_MAX):
        raise SystemExit("%s is %.2f x %.2f in; PLOS allows %.2f-%.2f in wide and at most %.2f in tall"
                         % (name, w_in, h_in, W_MIN, W_MAX, H_MAX))
    half = im.resize((im.size[0] * PNG_DPI // TIFF_DPI, im.size[1] * PNG_DPI // TIFF_DPI), Image.LANCZOS)
    return {"%s.tif" % name: encode(im, "tif", TIFF_DPI), "%s.png" % name: encode(half, "png", PNG_DPI)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare bytes; write nothing")
    ap.add_argument("--only", nargs="*", default=None, help="subset of figures to build")
    args = ap.parse_args()
    names = args.only or list(FIGURES)
    unknown = [n for n in names if n not in FIGURES]
    if unknown:
        raise SystemExit("unknown figures: %s" % unknown)
    OUT.mkdir(parents=True, exist_ok=True)
    bad = []
    for name in names:
        for fn, data in build(name).items():
            path = OUT / fn
            if args.check:
                if not path.exists() or path.read_bytes() != data:
                    bad.append(fn)
            else:
                path.write_bytes(data)
                print("wrote %s (%.2f MB)" % (path.relative_to(REPO), len(data) / 1e6))
    if args.check:
        if bad:
            print("FAIL: differ from a fresh build: %s" % ", ".join(bad))
            return 1
        print("PASS: %d figures reproduce byte for byte" % len(names))
        return 0
    if args.only is None:
        lines = []
        for p in sorted(OUT.glob("*")):
            if p.suffix in (".tif", ".png", ".json"):
                lines.append("%s  %s" % (hashlib.sha256(p.read_bytes()).hexdigest(), p.name))
        (OUT / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
