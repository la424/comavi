#!/usr/bin/env python3
"""The benchmark cohort flow figure (v10; superseded at v11 by
fig2_benchmark_detection.py and no longer built).

WHY THIS EXISTS
---------------
Readers consistently lost track of which subset a number referred to. The paper
names five populations -- 61 scored, 56 graded, 46 prioritized, 47 compared with
a pathogenicity predictor, 99 committed axis expectations -- and prose alone
cannot show that four of them are nested while the fifth is not.

That non-nesting is the actual source of the confusion and the thing this figure
has to make visible. The 47-variant pathogenicity set is NOT a subset of the 46
or the 56: it is the 49 complex variants minus the two with no AlphaMissense
score, so it KEEPS three variants that grading excludes and DROPS all twelve
single-chain variants. A reader following one chain of subsets hits that and the
chain breaks.

Every count is read from the committed records at build time and the branches are
asserted to sum, so the figure cannot drift from the table it describes.
"""
import pathlib

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO = pathlib.Path(__file__).resolve().parents[2]
RO = REPO / "reference_outputs"

GRADE = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
SILENT = "structurally_silent"


def cohort_numbers():
    """Every number in the figure, derived from the committed records."""
    can = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
    am = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_common_set.csv")
    can["graded"] = can.mech_consistency_t25.isin(GRADE)
    # The single-chain system is identified by having no partner, not by name.
    single = can.system == "brca1_brct"
    cx, br = can[~single], can[single]

    def split(d):
        g = d[d.graded]
        return dict(n=len(d), graded=len(g),
                    disruptive=int((g.expected_mech_class != SILENT).sum()),
                    nondisruptive=int((g.expected_mech_class == SILENT).sum()),
                    excluded=sorted(set(d.variant) - set(g.variant)))

    N = dict(total=len(can), systems=int(can.system.nunique()),
             cx=split(cx), br=split(br),
             cx_systems=int(cx.system.nunique()),
             graded=int(can.graded.sum()),
             axes=len(pd.read_csv(RO / "COMAVI_evidence_ledger.csv")),
             am_n=len(am), am_path=int((am.clinical_y == 1).sum()),
             am_benign=int((am.clinical_y == 0).sum()))
    N["disruptive"] = N["cx"]["disruptive"] + N["br"]["disruptive"]
    N["nondisruptive"] = N["cx"]["nondisruptive"] + N["br"]["nondisruptive"]
    # The non-nesting, shown as arithmetic rather than asserted in prose. The
    # comparison set is every variant carrying BOTH a clinical label and an
    # AlphaMissense score, which cuts across grading in two directions at once:
    # it keeps variants grading excluded, and drops graded ones whose score was
    # never retrieved. Stating that it "is not a subset" is weaker than showing
    # the two terms that make it one.
    g = set(can[can.graded].variant)
    av = set(am.variant)
    N["am_from_graded"] = len(av & g)
    N["am_from_excluded"] = len(av - g)
    N["am_dropped"] = sorted((set(cx.variant) & g) - av)
    N["am_dropped_single"] = len(set(br.variant) & av)
    assert N["am_from_graded"] + N["am_from_excluded"] == N["am_n"]
    assert N["cx"]["graded"] - len(N["am_dropped"]) == N["am_from_graded"]

    assert N["cx"]["n"] + N["br"]["n"] == N["total"]
    assert N["cx"]["graded"] + N["br"]["graded"] == N["graded"]
    assert N["disruptive"] + N["nondisruptive"] == N["graded"]
    assert N["am_path"] + N["am_benign"] == N["am_n"]
    assert N["am_dropped_single"] == 0, "the single-chain system has no clinical labels"

    # Internal arithmetic alone is not enough. Blanking a grade lowers both sides
    # of every sum above consistently, so the figure would draw a cohort the
    # paper does not describe and no assertion would fire. These tie the drawn
    # counts to the same committed records the Results text reads.
    import json
    sweep = {r["threshold_tag"]: r for r in
             json.loads((RO / "COMAVI_threshold_sweep.json").read_text())["rows"]}
    pv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
    assert N["graded"] == sweep["t25"]["whole_variant_n"], (
        "graded cohort %d disagrees with the threshold sweep denominator %d"
        % (N["graded"], sweep["t25"]["whole_variant_n"]))
    assert N["cx"]["graded"] == len(pv), (
        "prioritization track %d disagrees with the per-variant record %d"
        % (N["cx"]["graded"], len(pv)))
    assert N["cx"]["disruptive"] == int(pv.structural_ground_truth.astype(bool).sum())
    return N


# Geometry, in axes fractions. Five columns: the cohort splits into two parallel
# tracks which rejoin in a terminal box, so the figure opens and closes on a
# single number and the reader never has to add the tracks up themselves. A
# brace was tried first and drew through the boxes it was meant to group.
COL = [0.062, 0.232, 0.416, 0.608, 0.845]
WID = [0.118, 0.160, 0.156, 0.190, 0.150]
TOP, BOT, PATHO = 0.800, 0.545, 0.145
HEAD_Y = 0.962
BOXH, EXCH, PATH_H = 0.195, 0.165, 0.150


def _box(ax, cx, cy, w, h, lines, fc, ec, lw=0.7, ls="solid", zorder=2):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0,rounding_size=0.012",
        fc=fc, ec=ec, lw=lw, linestyle=ls, zorder=zorder))
    n = len(lines)
    step = 0.040
    y0 = cy + (n - 1) * step / 2
    for i, (txt, size, weight, colour) in enumerate(lines):
        ax.text(cx, y0 - i * step, txt, ha="center", va="center",
                fontsize=size, fontweight=weight, color=colour, zorder=zorder + 1)


def _arrow(ax, xy_from, xy_to, colour, ls="solid", lw=0.8):
    ax.add_patch(FancyArrowPatch(
        xy_from, xy_to, arrowstyle="-|>", mutation_scale=7, lw=lw,
        color=colour, linestyle=ls, shrinkA=0, shrinkB=0, zorder=1))


def fig_cohort(style, FULL, INK, GREY, RULE, STRUCT, NOLES, BAND):
    """Draw the cohort flow. Style constants are passed in so the figure cannot
    drift from the rest of the submission set."""
    N = cohort_numbers()
    style()
    fig = plt.figure(figsize=(FULL, 4.05))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Headers are left-aligned on each column's left edge rather than centred:
    # centring the first one ran it off the canvas, and left alignment also reads
    # as a column head rather than as a label belonging to the box below it.
    for cx, w, head in zip(COL[:4], WID[:4], ["Variants scored", "Structural context",
                                              "Not graded", "Graded by context"]):
        ax.text(cx - w / 2, HEAD_Y, head, ha="left", va="center", fontsize=8.5,
                fontweight="bold", color=GREY)

    mid = (TOP + BOT) / 2
    _box(ax, COL[0], mid, WID[0], BOXH,
         [("%d variants" % N["total"], 9.5, "bold", INK),
          ("%d protein" % N["systems"], 8, "normal", INK),
          ("systems", 8, "normal", INK)], "white", INK, lw=0.9)

    tracks = [
        (TOP, N["cx"], "%d variants" % N["cx"]["n"],
         ["%d multi-chain" % N["cx_systems"], "complexes"],
         ["no committed", "expectation"], "Prioritization set"),
        (BOT, N["br"], "%d variants" % N["br"]["n"],
         ["BRCA1 BRCT,", "single chain"],
         ["mechanism outside", "the model"], "Graded, not prioritized"),
    ]
    for y, d, title, sub, why, role in tracks:
        _box(ax, COL[1], y, WID[1], BOXH,
             [(title, 9.5, "bold", INK)] + [(s, 8, "normal", INK) for s in sub],
             "white", INK, lw=0.9)
        _box(ax, COL[2], y, WID[2], EXCH,
             [("%d excluded" % len(d["excluded"]), 8.5, "bold", GREY)]
             + [(s, 8, "normal", GREY) for s in why], BAND, RULE, lw=0.6)
        _box(ax, COL[3], y, WID[3], BOXH,
             [("%d graded" % d["graded"], 9.5, "bold", INK),
              ("%d disruptive" % d["disruptive"], 8, "normal", STRUCT),
              ("%d non-disruptive" % d["nondisruptive"], 8, "normal", NOLES)],
             "white", INK, lw=0.9)
        ax.text(COL[3], y - BOXH / 2 - 0.042, role, ha="center", va="center",
                fontsize=8, style="italic", color=GREY)
        _arrow(ax, (COL[0] + WID[0] / 2, mid), (COL[1] - WID[1] / 2, y), INK)
        _arrow(ax, (COL[1] + WID[1] / 2, y), (COL[2] - WID[2] / 2, y), INK)
        _arrow(ax, (COL[2] + WID[2] / 2, y), (COL[3] - WID[3] / 2, y), INK)
        _arrow(ax, (COL[3] + WID[3] / 2, y), (COL[4] - WID[4] / 2, mid), INK)

    _box(ax, COL[4], mid, WID[4], BOXH,
         [("%d graded" % N["graded"], 9.5, "bold", INK),
          ("%d disruptive" % N["disruptive"], 8, "normal", STRUCT),
          ("%d non-disruptive" % N["nondisruptive"], 8, "normal", NOLES)],
         "white", INK, lw=0.9)
    ax.text(COL[4], mid - BOXH / 2 - 0.042, "Mechanism-grading set",
            ha="center", va="center", fontsize=8, style="italic", color=GREY)

    # The pathogenicity comparison set, drawn below the rule because it is the
    # one population that is NOT a step in the flow above. Two feeds, one from
    # each side of the grading decision, are what make it non-nested, so it is
    # shown as arithmetic rather than asserted.
    ax.plot([0.012, 0.988], [0.312, 0.312], lw=0.6, ls=(0, (4, 3)), color=RULE)
    _box(ax, COL[4], PATHO, 0.190, PATH_H,
         [("%d variants" % N["am_n"], 9.5, "bold", INK),
          ("%d pathogenic," % N["am_path"], 8, "normal", INK),
          ("%d benign" % N["am_benign"], 8, "normal", INK)],
         "white", INK, lw=0.9, ls=(0, (3, 2)))
    ax.text(COL[4], PATHO - PATH_H / 2 - 0.042, "Pathogenicity comparison",
            ha="center", va="center", fontsize=8, style="italic", color=GREY)
    for src_y, label in [(0.212, "%d of the %d graded" % (N["am_from_graded"], N["cx"]["graded"])),
                         (0.078, "all %d excluded" % N["am_from_excluded"])]:
        _arrow(ax, (0.616, src_y), (COL[4] - 0.095, PATHO), GREY, ls=(0, (2.5, 2)), lw=0.7)
        ax.text(0.608, src_y, label, ha="right", va="center", fontsize=8, color=GREY)
    ax.text(0.012, PATHO,
            "Not a step in the flow above. It is\n"
            "every variant carrying both a clinical\n"
            "label and an AlphaMissense score, so it\n"
            "crosses the grading decision in both\n"
            "directions: it keeps the %d excluded,\n"
            "and drops %d graded variants with no\n"
            "score and all %d single-chain variants,\n"
            "which carry no clinical label."
            % (N["am_from_excluded"], len(N["am_dropped"]), N["br"]["n"]),
            ha="left", va="center", fontsize=8, color=GREY, linespacing=1.45)
    return fig, N
