#!/usr/bin/env python3
"""Supplementary figure: benchmark stress tests.

Reads the draws written by verification/stress_tests.py. Run that first:

    python3 verification/stress_tests.py \
        --canonical reference_outputs/scored_61var_canonical.csv \
        --scripts-dir scripts --out-dir verification_output
    python3 figures/src/figureS_stress_tests.py

Panel (c) is drawn as a discrete stem plot, not a histogram: mechanism-
consistency on n=57 graded variants moves in steps of 0.5/57 = 0.0088, so
binning invents gaps between attainable values.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FOCAL, COMP, GREY = "#B02A3A", "#1B7A5A", "#6E6E6E"


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "figure.dpi": 150,
        "savefig.bbox": None, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", default="verification_output")
    ap.add_argument("--out-dir", default="figures")
    args = ap.parse_args()

    ind, outd = Path(args.in_dir), Path(args.out_dir)
    z = np.load(ind / "comavi_stress_draws.npz")
    loso = pd.read_csv(ind / "comavi_leave_one_system_out.csv").sort_values("mc")

    mc_null = z["mc_null"]
    mc_obs = float(z["mc_obs"][0])
    p_mc = float(z["p_mc"][0])
    noise_mc, noise_flips = z["noise_mc"], z["noise_flips"]

    style()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))

    # (a) permutation null
    ax = axes[0]
    ax.hist(mc_null, bins=28, color="#BFBFBF", edgecolor="none")
    ymax = ax.get_ylim()[1]
    ax.axvline(mc_obs, color=FOCAL, lw=1.8)
    ax.annotate(f"observed {mc_obs:.2f}", xy=(mc_obs, ymax * 0.50),
                xytext=(mc_obs - 0.050, ymax * 0.63), color=FOCAL, ha="right",
                fontsize=6, arrowprops=dict(arrowstyle="-", color=FOCAL, lw=0.7))
    ax.text(np.percentile(mc_null, 15), ymax * 0.86, "shuffled\nground truth",
            ha="center", fontsize=6, color=GREY, linespacing=1.2)
    exp = int(np.floor(np.log10(p_mc)))
    ax.text(mc_obs - 0.11, ymax * 0.97,
            rf"$p$ = {p_mc / 10**exp:.0f}$\times$10$^{{{exp}}}$",
            ha="center", va="top", fontsize=6, color=FOCAL)
    ax.set_xlabel("Mechanism-consistency")
    ax.set_ylabel("Permutations")
    ax.set_title("Far above chance", loc="left")
    ax.margins(x=0.04)

    # (b) leave-one-system-out
    ax = axes[1]
    y = np.arange(len(loso))
    ax.hlines(y, mc_obs, loso.mc, color="#CCCCCC", lw=0.8, zorder=1)
    ax.scatter(loso.mc, y, s=14, color=FOCAL, zorder=3)
    ax.axvline(mc_obs, color="#444444", lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.text(mc_obs - 0.0018, -1.45, f"full set {mc_obs:.2f}", fontsize=6,
            color="#444444", va="bottom", ha="right")
    ax.set_yticks(y)
    ax.set_yticklabels([g.upper() for g in loso.dropped], fontsize=6)
    span = loso.mc.max() - loso.mc.min()
    ax.set_xlim(loso.mc.min() - 0.18 * span, loso.mc.max() + 0.08 * span)
    ax.set_ylim(-1.6, len(loso) - 0.4)
    ax.set_xlabel("Mechanism-consistency, one system dropped")
    ax.set_title("No system carries the result", loc="left")

    # (c) replicate noise, discrete
    ax = axes[2]
    vals, cnts = np.unique(np.round(noise_mc, 6), return_counts=True)
    ax.vlines(vals, 0, cnts, color=COMP, lw=3.2, zorder=2)
    ax.scatter(vals, cnts, s=9, color=COMP, zorder=3)
    ax.axvline(mc_obs, color=FOCAL, lw=1.4, zorder=1)
    ax.text(mc_obs - 0.0030, cnts.max() * 0.62, "reported", fontsize=6,
            color=FOCAL, ha="center", va="center", rotation=90)
    ax.text(0.98, 0.97,
            f"sd = {noise_mc.std():.3f}\n"
            f"{100 * np.mean(noise_flips == 0):.0f}% of draws\nlabel-identical",
            transform=ax.transAxes, ha="right", va="top", fontsize=6,
            color=COMP, linespacing=1.25)
    # Pad to the next tick beyond the data so the rightmost tick label is not
    # flush with the panel edge (it overhangs the canvas at right=0.995).
    ax.set_xlim(vals.min() - 0.010, vals.max() + 0.040)
    ax.set_xlabel("Mechanism-consistency")
    ax.set_ylabel("Perturbation draws")
    ax.set_title("Stable to force-field noise", loc="left")
    ax.margins(y=0.08)

    # v7.3: panel (b) y-tick labels are full system names and extend left into
    # the inter-panel gap; pooling also moved the observed MC to the right edge
    # of panel (a), so its annotation now reaches the same gap. Widen it.
    fig.subplots_adjust(left=0.075, right=0.995, top=0.80, bottom=0.20, wspace=0.58)
    for a, lab in zip(axes, "abc"):
        bb = a.get_position()
        fig.text(bb.x0 - 0.060, bb.y1 + 0.075, lab, fontsize=9,
                 fontweight="bold", va="top", ha="left")

    # Geometric self-check: no label may overlap another label, cross a spine,
    # or fall outside the canvas. Fails loudly rather than shipping silently.
    # Collect labels explicitly: a fig.findobj(Text) walk also returns the spare
    # tick artists matplotlib keeps outside the view limits, which produce
    # phantom overlaps.
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    def visible_ticks(axis):
        """Tick labels actually inside the view interval."""
        lo, hi = sorted(axis.get_view_interval())
        return [t for t, loc in zip(axis.get_ticklabels(), axis.get_ticklocs())
                if lo - 1e-9 <= loc <= hi + 1e-9
                and t.get_text().strip() and t.get_visible()]

    ticks, items = {}, list(fig.texts)
    for a in fig.axes:
        tk = visible_ticks(a.xaxis) + visible_ticks(a.yaxis)
        ticks[a] = set(tk)
        items += tk + list(a.texts) + [a.title, a.xaxis.label, a.yaxis.label]
    texts = [(t, t.get_window_extent(rend))
             for t in items if t.get_text().strip() and t.get_visible()]
    bad = [(u.get_text(), v.get_text())
           for i, (u, bu) in enumerate(texts) for v, bv in texts[i + 1:]
           if bu.overlaps(bv)]
    bad += [(t.get_text(), "spine") for t, bt in texts for a in fig.axes
            for s in a.spines.values()
            if s.get_visible() and bt.overlaps(s.get_window_extent(rend))
            and t not in ticks[a]]
    oob = [t.get_text() for t, b in texts
           if not (b.x0 >= -1 and b.y0 >= -1
                   and b.x1 <= fig.bbox.x1 + 1 and b.y1 <= fig.bbox.y1 + 1)]
    if bad or oob:
        raise SystemExit(f"FAIL layout check: overlaps={bad} out_of_bounds={oob}")
    print("layout check: no overlaps, nothing clipped")

    outd.mkdir(parents=True, exist_ok=True)
    fig.savefig(outd / "COMAVI_FigS_stress_tests.png", dpi=300)
    fig.savefig(outd / "COMAVI_FigS_stress_tests.pdf")
    print(f"MC observed {mc_obs:.4f} | null {mc_null.mean():.4f} | p {p_mc:.5f}")
    print(f"LOSO MC range {loso.mc.min():.4f}-{loso.mc.max():.4f}")
    print(f"noise sd {noise_mc.std():.4f} | "
          f"{100 * np.mean(noise_flips == 0):.0f}% label-identical")
    print(f"wrote {outd}/COMAVI_FigS_stress_tests.png/.pdf")


if __name__ == "__main__":
    main()
