#!/usr/bin/env python3
"""Fig 2: the benchmark populations, and what scoring the assembled complex adds.

WHY THIS EXISTS
---------------
v11 rebuilt Fig 2 around two reviewer comments. The populations comment asked
for one cohort that every claim nests inside: the benchmark is now the 56
graded variants (the five candidates no gradeable expectation could be defined
for are named once, in Methods), and with clinical labels classified by the v11 rule (ClinVar expert panel where one exists, otherwise Franklin)
the pathogenicity comparison is the labelled part of the prioritization
population, so the populations nest (56 > 46 > 39). Panel A draws that.

The detection comment asked for a figure in place of a paragraph of counts.
Panels B-D carry them:
  B  detection by curated class at 2.5 kcal/mol, split by the first axis that
     fires (isolated subunit, then assembled-complex stability, then binding)
     in the Fig 1 axis colours;
  C  the same split for detections and false calls at all five thresholds;
  D  detection, localization and false calls when axes are removed from the
     analysis (COMAVI's own energies, recounted; a comparison of designs, not
     of the published tools).

Every count is read from the committed records and the panels are asserted to
agree with each other and with the threshold sweep, so the figure cannot drift
from the text that cites it.
"""
import json
import pathlib
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

REPO = pathlib.Path(__file__).resolve().parents[2]
RO = REPO / "reference_outputs"
sys.path.insert(0, str(REPO / "scripts"))
import benchmark_population as bp  # noqa: E402

SILENT = "structurally_silent"
SINGLE_CHAIN = "brca1_brct"
CLASSES = [("Fold", ("fold_mechanism",)),
           ("Interface", ("ppi_destab_mechanism", "ppi_stab_mechanism")),
           ("Mixed", ("mixed_structural",))]
AXES = ("isolated_subunit", "assembled_stability", "binding_only")
TAGS = ["t10", "t15", "t20", "t25", "tSAP"]


def J(name):
    return json.loads((RO / name).read_text())


def numbers():
    """Every count in the figure, from the committed records, cross-checked."""
    can = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
    b = bp.benchmark(can)
    silent = b.expected_mech_class == SILENT
    single = b.system == SINGLE_CHAIN
    cx, br = b[~single], b[single]
    pv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
    am = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_common_set.csv")
    dva = J("COMAVI_detection_vs_attribution.json")
    sweep = {r["threshold_tag"]: r for r in J("COMAVI_threshold_sweep.json")["rows"]}

    N = dict(n=len(b), systems=int(b.system.nunique()),
             mech=int((~silent).sum()), nondis=int(silent.sum()),
             cx=len(cx), cx_systems=int(cx.system.nunique()),
             cx_mech=int((cx.expected_mech_class != SILENT).sum()),
             cx_nondis=int((cx.expected_mech_class == SILENT).sum()),
             br=len(br), br_mech=int((br.expected_mech_class != SILENT).sum()),
             br_nondis=int((br.expected_mech_class == SILENT).sum()),
             am=len(am), am_path=int((am.clinical_y == 1).sum()),
             am_benign=int((am.clinical_y == 0).sum()))
    cx_keys = set(zip(cx.system, cx.variant))
    am_keys = set(zip(am.system, am.variant))
    not_compared = cx[[k not in am_keys for k in zip(cx.system, cx.variant)]]
    N["not_compared"] = len(not_compared)
    N["not_compared_labels"] = not_compared.phenotype.value_counts().to_dict()

    # The nesting the panel draws, and the records the text reads.
    assert N["n"] == bp.N_BENCHMARK == sweep["t25"]["whole_variant_n"]
    assert N["cx"] + N["br"] == N["n"] and N["cx_mech"] + N["br_mech"] == N["mech"]
    assert N["cx"] == len(pv) and N["cx_mech"] == int(pv.structural_ground_truth.astype(bool).sum())
    assert am_keys <= cx_keys, "the pathogenicity comparison must nest in the prioritization population"
    assert N["am"] + N["not_compared"] == N["cx"]
    assert set(N["not_compared_labels"]) <= {"vus", "unclassified"}, N["not_compared_labels"]
    assert N["mech"] == dva["detection"]["n"] and N["nondis"] == dva["correct_rejection"]["n"]

    pc = dva["per_class"]
    cls = []
    for label, members in CLASSES:
        n = sum(pc[m]["n"] for m in members if m in pc)
        split = {a: sum(pc[m]["first_firing_axis"][a] for m in members if m in pc) for a in AXES}
        assert sum(split.values()) == sum(pc[m]["detection"] for m in members if m in pc)
        cls.append((label, n, split))
    nd_split = pc[SILENT]["first_firing_axis"]
    assert sum(n for _, n, _ in cls) == N["mech"]
    assert sum(sum(s.values()) for _, _, s in cls) == dva["detection"]["k"]
    assert sum(nd_split.values()) == N["nondis"] - dva["correct_rejection"]["k"]
    assert sweep["t25"]["detection_first_axis"] == {a: sum(s[a] for _, _, s in cls) for a in AXES}
    N["classes"], N["nd_split"] = cls, nd_split
    N["sweep"] = [(sweep[t]["detection_first_axis"], sweep[t]["false_calls_first_axis"]) for t in TAGS]
    N["designs"] = dva["axis_designs"]
    assert N["designs"]["all_three"]["detection"] == dva["detection"]["k"]
    return N


def _box(ax, cx, cy, w, lines, ink, rule, step_in, lw=0.8):
    """A rounded box sized to its lines; `step_in` is the line pitch in inches."""
    h_ax = ax.get_position().height * ax.figure.get_figheight()
    step = step_in / h_ax
    h = step * (len(lines) + 0.55)
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=0.015",
                                fc="white", ec=ink, lw=lw, zorder=2))
    y0 = cy + (len(lines) - 1) * step / 2
    for i, (txt, size, weight, colour) in enumerate(lines):
        ax.text(cx, y0 - i * step, txt, ha="center", va="center", fontsize=size,
                fontweight=weight, color=colour, zorder=3)
    return h


def _arrow(ax, a, b, colour, lw=0.8, ls="solid"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=7, lw=lw,
                                 color=colour, linestyle=ls, shrinkA=0, shrinkB=0, zorder=1))


def panel_a(ax, N, C):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    step = 0.15
    X0, X1, X2 = 0.115, 0.455, 0.815
    W0, W1, W2 = 0.20, 0.235, 0.215
    TOP, BOT, MID = 0.735, 0.255, 0.495
    h0 = _box(ax, X0, MID, W0, [("%d variants" % N["n"], 9.5, "bold", C["INK"]),
                                 ("%d protein systems" % N["systems"], 8, "normal", C["INK"]),
                                 ("%d mechanisms" % N["mech"], 8, "normal", C["STRUCT"]),
                                 ("%d non-disruptive" % N["nondis"], 8, "normal", C["NOLES"])],
              C["INK"], C["RULE"], step, lw=0.9)
    h1 = _box(ax, X1, TOP, W1, [("%d variants" % N["cx"], 9.5, "bold", C["INK"]),
                                 ("%d multi-chain complexes" % N["cx_systems"], 8, "normal", C["INK"]),
                                 ("%d mechanisms" % N["cx_mech"], 8, "normal", C["STRUCT"]),
                                 ("%d non-disruptive" % N["cx_nondis"], 8, "normal", C["NOLES"])],
              C["INK"], C["RULE"], step)
    h2 = _box(ax, X1, BOT, W1, [("%d variants" % N["br"], 9.5, "bold", C["INK"]),
                                 ("BRCA1 BRCT, single chain", 8, "normal", C["INK"]),
                                 ("%d mechanisms" % N["br_mech"], 8, "normal", C["STRUCT"]),
                                 ("%d non-disruptive" % N["br_nondis"], 8, "normal", C["NOLES"])],
              C["INK"], C["RULE"], step)
    h3 = _box(ax, X2, TOP, W2, [("%d variants" % N["am"], 9.5, "bold", C["INK"]),
                                 ("Clinically classified", 8, "normal", C["INK"]),
                                 ("%d pathogenic" % N["am_path"], 8, "normal", C["INK"]),
                                 ("%d benign" % N["am_benign"], 8, "normal", C["INK"])],
              C["INK"], C["RULE"], step)
    for x, y, h, txt in [(X0, MID, h0, "Benchmark: all graded"),
                         (X1, TOP, h1, "Prioritization population"),
                         (X1, BOT, h2, "Graded; no partner, no tier"),
                         (X2, TOP, h3, "Pathogenicity comparison")]:
        ax.text(x, y - h / 2 - 0.055, txt, ha="center", va="center", fontsize=8,
                style="italic", color=C["GREY"])
    _arrow(ax, (X0 + W0 / 2, MID), (X1 - W1 / 2, TOP), C["INK"])
    _arrow(ax, (X0 + W0 / 2, MID), (X1 - W1 / 2, BOT), C["INK"])
    _arrow(ax, (X1 + W1 / 2, TOP), (X2 - W2 / 2, TOP), C["INK"])
    lab = N["not_compared_labels"]
    ax.text(X2 - W2 / 2, BOT + 0.02,
            "%d variants in complexes without a\npathogenic or benign classification\n"
            "(%d uncertain significance, %d unclassified)\nare not compared"
            % (N["not_compared"], lab.get("vus", 0), lab.get("unclassified", 0)),
            ha="left", va="center", fontsize=8, color=C["GREY"], linespacing=1.4)


def _legend_axes(fig, ax, handles, ncol, dy_in):
    tr = ax.transAxes + mtrans.ScaledTranslation(0, -dy_in, fig.dpi_scale_trans)
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, 0), bbox_transform=tr,
               ncol=ncol, frameon=False, handlelength=1.2, columnspacing=1.0,
               handletextpad=0.4, fontsize=8)


def panel_b(ax, N, C):
    rows = [(lab, n, split, False) for lab, n, split in N["classes"]]
    rows.append(("Non-disruptive", N["nondis"], N["nd_split"], True))
    ys = [3.6, 2.6, 1.6, 0.0]
    colours = [C["MONO"], C["FOLD"], C["BIND"]]
    for (lab, n, split, nd), y in zip(rows, ys):
        left = 0
        for a, col in zip(AXES, colours):
            if split[a]:
                ax.barh(y, split[a], left=left, height=0.62, color=col, lw=0)
            left += split[a]
        ax.barh(y, n - left, left=left, height=0.62, color=C["BAND"], lw=0,
                hatch="////" if nd else None, edgecolor="white" if nd else None)
        ax.text(n + 0.6, y, ("%d false calls" % left) if nd else ("%d of %d" % (left, n)),
                va="center", ha="left", fontsize=8, color=C["INK"])
    ax.set_yticks(ys, ["%s (%d)" % (lab, n) for lab, n, _, _ in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 39)
    ax.set_ylim(-0.6, 4.2)
    ax.set_xticks([0, 10, 20])
    ax.set_xlabel("Variants at 2.5 kcal/mol")
    ax.spines["left"].set_visible(False)


def panel_c(ax, N, C, xlab):
    colours = [C["MONO"], C["FOLD"], C["BIND"]]
    for i, (det, fc) in enumerate(N["sweep"]):
        up = down = 0
        for a, col in zip(AXES, colours):
            if det[a]:
                ax.bar(i, det[a], bottom=up, width=0.62, color=col, lw=0)
            up += det[a]
            if fc[a]:
                ax.bar(i, -fc[a], bottom=-down, width=0.62, color=col, lw=0, alpha=0.55)
            down += fc[a]
        ax.text(i, up + 0.6, "%d" % up, ha="center", va="bottom", fontsize=8, color=C["INK"])
        ax.text(i, -down - 0.6, "%d" % down, ha="center", va="top", fontsize=8, color=C["INK"])
    ax.axhline(0, color=C["INK"], lw=0.6)
    # The per-axis limit reads 2.9/2.9/3.5 everywhere else in the paper; it is
    # broken over two lines here only because the panel is narrow.
    ax.set_xticks(range(5), list(xlab[:-1]) + [xlab[-1].replace("/3.5", "/\n3.5")])
    ax.set_xlim(-0.55, 4.55)
    ax.set_ylim(-21, 31)
    ax.set_yticks([-15, -10, -5, 0, 5, 10, 15, 20, 25], ["15", "10", "5", "0", "5", "10", "15", "20", "25"])
    ax.set_xlabel("Energy threshold (kcal/mol)")
    ax.text(-0.27, 0.78, "Detected (of %d)" % N["mech"], transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=8.5)
    ax.text(-0.27, 0.20, "False calls (of %d)" % N["nondis"], transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=8.5)


def panel_d(ax, N, C):
    order = [("isolated_subunit", "Isolated\nsubunit"),
             ("isolated_subunit_and_binding", "Isolated +\nbinding"),
             ("isolated_subunit_and_assembled_stability", "Isolated +\nassembled\nstability"),
             ("all_three", "All three\n(COMAVI)")]
    bars = [("detection", C["INK"], None, "detected (of %d)" % N["mech"]),
            ("attribution", "#A6A6A6", None, "on the expected axis (of %d)" % N["mech"]),
            ("false_calls", "white", "////", "false calls (of %d)" % N["nondis"])]
    for j, (key, _lab) in enumerate(order):
        d = N["designs"][key]
        for k, (field, col, hatch, _) in enumerate(bars):
            y = (3 - j) * 1.0 + (1 - k) * 0.27
            ax.barh(y, d[field], height=0.24, color=col, hatch=hatch,
                    edgecolor=C["INK"] if hatch else col, lw=0.5 if hatch else 0)
            ax.text(d[field] + 0.5, y, "%d" % d[field], va="center", ha="left", fontsize=8)
    ax.set_yticks([3, 2, 1, 0], [lab for _, lab in order])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 25)
    ax.set_ylim(-0.55, 3.55)
    ax.set_xticks([0, 10, 20])
    ax.set_xlabel("Variants at 2.5 kcal/mol")
    ax.spines["left"].set_visible(False)
    return [Patch(fc=col, ec=C["INK"] if hatch else col, hatch=hatch, lw=0.5 if hatch else 0,
                  label=lab) for _, col, hatch, lab in bars]


def fig2(style, letter, FULL, C, xlab):
    N = numbers()
    style()
    fig = plt.figure(figsize=(FULL, 5.35))
    ax_a = fig.add_axes([0.0, 0.615, 1.0, 0.385])
    ax_b = fig.add_axes([0.165, 0.175, 0.20, 0.33])
    ax_c = fig.add_axes([0.49, 0.175, 0.235, 0.33])
    ax_d = fig.add_axes([0.83, 0.175, 0.155, 0.33])
    panel_a(ax_a, N, C)
    panel_b(ax_b, N, C)
    panel_c(ax_c, N, C, xlab)
    d_handles = panel_d(ax_d, N, C)
    for ax, s, dx in [(ax_a, "A", 0.0), (ax_b, "B", -1.0), (ax_c, "C", -0.45), (ax_d, "D", -0.62)]:
        letter(ax, s, dx=dx, dy=0.06 if ax is not ax_a else -0.2)
    axis_handles = [Patch(fc=C["MONO"], label="Isolated subunit"),
                    Patch(fc=C["FOLD"], label="Assembled-complex stability"),
                    Patch(fc=C["BIND"], label="Binding only"),
                    Patch(fc=C["BAND"], label="Not reached"),
                    Patch(fc=C["BAND"], hatch="////", ec="white", label="Correctly rejected")]
    _legend_axes(fig, ax_b, axis_handles, 3, 0.43)
    _legend_axes(fig, ax_d, d_handles, 1, 0.43)
    return fig, N
