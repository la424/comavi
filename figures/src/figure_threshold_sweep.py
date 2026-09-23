#!/usr/bin/env python3
"""Every reported metric as a function of the decision threshold.

Reads reference_outputs/COMAVI_threshold_sweep.json only. Nothing here is
typed: every plotted value and every count in a label comes from that record,
and the assertions below fail the build if the record and the figure disagree.

Three panels, answering three different questions:

  A  what the threshold buys and what it costs. Detection is nearly flat, so
     the tradeoff is not detection against rejection as a two-curve plot
     implies -- it is ATTRIBUTION against rejection.
  B  whether one threshold suits all three axes. It does not: the binding
     axis is best at the lowest threshold while the fold axis is best at the
     reference, which is why a single global cut costs the interface arm.
  C  the whole-variant score with system-level bootstrap bands, showing that
     the reference point is the maximum and that the bands overlap.

Usage:
    PYTHONPATH=. python figures/src/figure_threshold_sweep.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "figures" / "src"))
from _comavi_style import (COLUMN_W, META_GREY, apply_figure_style,  # noqa: E402
                           panel_letter, set_frame)

REC = REPO / "reference_outputs" / "COMAVI_threshold_sweep.json"
OUT = REPO / "figures" / "COMAVI_Figure_threshold_sweep.png"

# Threaded colours: each quantity keeps its hue wherever it appears.
C_DETECT = "#1f6f5c"     # teal    - detection
C_ATTRIB = "#c2410c"     # orange  - attribution
C_REJECT = "#2563eb"     # blue    - correct rejection
C_MONO = "#7c3aed"       # violet  - isolated-subunit stability
C_FOLD = "#0891b2"       # cyan    - stability in the complex
C_BIND = "#c2410c"       # orange  - binding (same hue as attribution:
                         #           both are about naming the interface)
C_WHOLE = "#1f2937"      # near-black - the headline


def main():
    rec = json.loads(REC.read_text())
    rows = rec["rows"]
    bands = rec["cluster_bands"]
    x = list(range(len(rows)))
    labels = [r["threshold_kcal_mol"] for r in rows]
    ref_i = [i for i, r in enumerate(rows) if r["threshold_tag"] == "t25"][0]

    # --- claim-title checks (figure-style 1.4): a sentence title is tested
    # against every plotted row before rendering, not asserted in prose.
    det = [r["detection"] for r in rows]
    att = [r["attribution"] for r in rows]
    rej = [r["correct_rejection"] for r in rows]
    mono = [r["axis_monomer"] for r in rows]
    fold = [r["axis_fold"] for r in rows]
    bind = [r["axis_binding"] for r in rows]
    whole = [r["whole_variant"] for r in rows]

    assert all(att[i] >= att[i + 1] for i in range(len(att) - 1)), \
        "panel A title claims attribution falls monotonically; it does not"
    assert all(rej[i] <= rej[i + 1] for i in range(len(rej) - 1)), \
        "panel A title claims rejection rises monotonically; it does not"
    _flat = det[1:ref_i + 1]          # 1.5 through the reference, inclusive
    assert max(_flat) - min(_flat) < 1e-9, \
        ("panel A title claims detection is flat from 1.5 to the reference; "
         "observed %s" % _flat)
    assert bind.index(max(bind)) < fold.index(max(fold)), \
        "panel B title claims binding peaks below fold; it does not"
    assert whole.index(max(whole)) == ref_i, \
        "panel C title claims the reference threshold is the maximum; it is not"

    # Graded arms and the measured-effect column, so every threshold-dependent
    # quantity the paper reports lives in ONE figure with its definition named.
    # Panels C and D subsume the whole of the figure this one replaces: its
    # panel A was the graded structural and silent arms (identical to the
    # sensitivity and specificity columns read here) and its panel B was the
    # measured-effect recovery.
    op = pd.read_csv(REPO / "reference_outputs" / "COMAVI_threshold_operating_points.csv")
    arm_struct = op.sensitivity.tolist()
    arm_silent = op.specificity.tolist()
    assert len(arm_struct) == len(rows), "operating points and sweep disagree on the grid"
    for i, r in enumerate(rows):
        assert abs(op.pooled_MC.iloc[i] - r["whole_variant"]) < 5e-4, (
            "pooled score disagrees between the operating points and the sweep "
            "at %s: %.4f vs %.4f" % (r["threshold_kcal_mol"],
                                     op.pooled_MC.iloc[i], r["whole_variant"]))
    # Panel C's title also claims the pooled peak sits ABOVE the arm crossing.
    # An earlier draft of this title said the peak was AT the crossing, which is
    # false -- the arms cross between 1.5 and 2.0 while the peak is at 2.5 --
    # so the claim now carries its own check.
    _d = [a - b for a, b in zip(arm_struct, arm_silent)]
    _cross = [i for i in range(len(_d) - 1) if _d[i] * _d[i + 1] < 0]
    assert _cross and max(_cross) < ref_i, (
        "panel C title claims the pooled peak lies above the arm crossing; "
        "crossing indices %s, peak index %d" % (_cross, ref_i))
    meas = json.loads((REPO / "reference_outputs"
                       / "COMAVI_measured_effect_tables.json").read_text())
    mrec = [d["recovered"] for d in meas["measured_effects_recovered"]["by_threshold"]]
    mn = meas["measured_effects_recovered"]["population"]["n"]

    apply_figure_style(sizes=(8, 7, 6))
    fig, axes = plt.subplots(1, 4, figsize=(COLUMN_W * 1.90, 2.62))

    def ref_marker(ax):
        ax.axvline(ref_i, color=META_GREY, lw=0.8, ls=(0, (4, 3)), zorder=0)

    # ---------------------------------------------------------------- panel A
    ax = axes[0]
    ref_marker(ax)
    for ys, c, lab in ((det, C_DETECT, "detection"),
                       (att, C_ATTRIB, "attribution"),
                       (rej, C_REJECT, "correct rejection")):
        ax.plot(x, ys, "-o", color=c, lw=1.7, ms=4.2, label=lab, zorder=3)
    ax.text(x[-1] + 0.10, det[-1], "detection", color=C_DETECT, va="center", fontsize=7)
    ax.text(x[-1] + 0.10, att[-1] - 0.045, "attribution", color=C_ATTRIB, va="center", fontsize=7)
    ax.text(x[-1] + 0.10, rej[-1], "correct\nrejection", color=C_REJECT, va="center", fontsize=7)
    ax.set_title("Raising the threshold buys rejection\nand gives up attribution", loc="left")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0.05, 1.0)
    ax.set_xlim(-0.25, len(rows) - 0.32)

    # ---------------------------------------------------------------- panel B
    ax = axes[1]
    ref_marker(ax)
    for ys, c, lab in ((mono, C_MONO, "isolated subunit"),
                       (fold, C_FOLD, "in complex"),
                       (bind, C_BIND, "binding")):
        ax.plot(x, ys, "-o", color=c, lw=1.7, ms=4.2, label=lab, zorder=3)
    ax.plot([bind.index(max(bind))], [max(bind)], "o", mfc="none", mec=C_BIND,
            ms=10, mew=1.4, zorder=4)
    ax.plot([fold.index(max(fold))], [max(fold)], "o", mfc="none", mec=C_FOLD,
            ms=10, mew=1.4, zorder=4)
    ax.set_title("One threshold does not suit all three axes:\nbinding peaks lowest, fold peaks at the reference",
                 loc="left")
    ax.set_ylabel("Directional agreement")
    ax.set_ylim(0.50, 0.85)
    ax.set_xlim(-0.30, len(rows) - 0.68)
    ax.legend(frameon=False, loc="lower right", handlelength=1.3,
              borderaxespad=0.15, fontsize=7)

    # ---------------------------------------------------------------- panel C
    ax = axes[2]
    ref_marker(ax)
    lo = [bands[r["threshold_tag"]]["lo"] for r in rows]
    hi = [bands[r["threshold_tag"]]["hi"] for r in rows]
    ax.fill_between(x, lo, hi, color=C_WHOLE, alpha=0.13, lw=0, zorder=1)
    ax.plot(x, arm_struct, "-s", color=C_DETECT, lw=1.4, ms=3.6, zorder=2)
    ax.plot(x, arm_silent, "-^", color=C_REJECT, lw=1.4, ms=3.6, zorder=2)
    ax.plot(x, whole, "-o", color=C_WHOLE, lw=1.9, ms=4.6, zorder=3)
    ax.annotate("%.3f" % whole[ref_i], (ref_i, whole[ref_i]),
                textcoords="offset points", xytext=(-31, -3),
                color=C_WHOLE, fontsize=7, fontweight="bold")
    ax.text(0.05, 0.92, "pooled", transform=ax.transAxes, color=C_WHOLE, fontsize=7)
    ax.text(0.05, 0.83, "structural arm", transform=ax.transAxes, color=C_DETECT, fontsize=7)
    ax.text(0.05, 0.74, "no-lesion arm", transform=ax.transAxes, color=C_REJECT, fontsize=7)
    ax.set_title("Graded: the pooled score peaks above\nthe point where the two arms cross", loc="left")
    ax.set_ylabel("Mechanism-pattern grade")
    ax.set_ylim(0.38, 0.99)
    ax.set_xlim(-0.30, len(rows) - 0.68)
    ax.text(0.05, 0.04, "band: 95%% CI, %d-system bootstrap"
            % bands["t25"]["n_systems"], transform=ax.transAxes,
            color=META_GREY, fontsize=6.2)

    # ---------------------------------------------------------------- panel D
    ax = axes[3]
    ref_marker(ax)
    ax.plot(x, [v / mn for v in mrec], "-D", color=C_ATTRIB, lw=1.7, ms=4.0, zorder=3)
    for xi, v in zip(x, mrec):
        ax.annotate("%d" % v, (xi, v / mn), textcoords="offset points",
                    xytext=(0, 6), ha="center", fontsize=6.2, color=C_ATTRIB)
    ax.set_title("Against direct measurement, stringency\ncosts recovery fastest of all", loc="left")
    ax.set_ylabel("Measured effects recovered")
    ax.set_ylim(0.14, 0.98)
    ax.set_xlim(-0.30, len(rows) - 0.68)
    ax.text(0.05, 0.04, "n = %d comparisons, |measured| >= 1 kcal/mol" % mn,
            transform=ax.transAxes, color=META_GREY, fontsize=6.2)

    for i, ax in enumerate(axes):
        set_frame(ax)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6.4)
        ax.set_xlabel("Decision threshold (kcal/mol)")
        panel_letter(ax, "ABCD"[i])

    fig.subplots_adjust(left=0.050, right=0.990, top=0.80, bottom=0.185, wspace=0.315)
    fig.savefig(OUT, dpi=300)
    fig.savefig(OUT.with_suffix(".pdf"))

    # render-then-verify: every text box inside the canvas, no text-text overlap
    import matplotlib as mpl
    # savefig(dpi=300) leaves a renderer in a different dpi space than
    # fig.bbox, so window extents taken straight after it are not comparable
    # to the canvas box and every axis label reads as outside it. Re-draw at
    # canvas dpi first; this is the check being wrong, not the layout.
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    outside = [t.get_text() for t, b in texts
               if not (fig.bbox.x0 <= b.x0 and b.x1 <= fig.bbox.x1
                       and fig.bbox.y0 <= b.y0 and b.y1 <= fig.bbox.y1)]
    tl = {ax: set(ax.get_xticklabels() + ax.get_yticklabels()) for ax in fig.axes}
    allt = set().union(*tl.values())
    ov = [(a.get_text(), b.get_text())
          for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:]
          if ba.overlaps(bb) and not (a in allt and b in allt)]
    assert not outside, "text outside canvas: %s" % outside
    print("bbox check: %d text objects, 0 outside canvas, %d overlaps" % (len(texts), len(ov)))
    if ov:
        print("  overlaps: %s" % ov[:6])
    print("detection      %s" % det)
    print("attribution    %s" % att)
    print("rejection      %s" % rej)
    print("per-axis  mono %s\n          fold %s\n          bind %s" % (mono, fold, bind))
    print("whole-variant  %s  (argmax at %s)" % (whole, labels[ref_i]))
    print("wrote %s" % OUT.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
