"""Shared canvas geometry and font ladder for COMAVI manuscript figures.

Why this module exists
----------------------
The figures were originally authored at whatever canvas width suited each
one (7.2, 9.2, 11.6, 13.4, 13.8 in). Every figure then lands in the same
6.5-inch manuscript text column, so each was scaled by a *different* factor
(0.90x down to 0.47x) and its type shrank by that same factor. Measured on
the v22 embeds, effective type ran from 7.3 pt (fine) down to 3.7 pt
(illegible, and below any journal's floor).

The fix is not to enlarge fonts -- it is to author every figure at ONE
width, so that the page scale factor is the same for all of them and the
authored point size is the delivered point size.

    authored at COLUMN_W = 7.2 in  ->  placed at DOC_W = 6.5 in
    scale 0.903, so an 8.0 pt label is delivered at 7.2 pt.

Font ladder
-----------
figure-style Sec 5.2: at most three sizes, mapped to ROLE not to available
space, plus the panel letter as the single permitted exception. If a label
does not fit at its role's size, the layout gets fixed -- not the size.

    BASE   8.0  titles, axis labels, series identity / direct labels
    SECOND 7.0  legend text, annotations, on-mark values
    TICK   6.0  tick labels
    LETTER 9.0  panel letters (bold)

Delivered at 0.903x these are 7.2 / 6.3 / 5.4 pt. The 5.4 pt tick floor is
the reason no figure may be authored wider than COLUMN_W.
"""
import matplotlib as mpl

# ---- geometry -------------------------------------------------------------
COLUMN_W = 7.2   # authored width for every manuscript figure, inches
DOC_W = 6.5      # width of the Word/manuscript text column, inches
PAGE_SCALE = DOC_W / COLUMN_W   # 0.903

# ---- font ladder ----------------------------------------------------------
BASE, SECOND, TICK = 8.0, 7.0, 6.0
LETTER = 9.0

META_GREY = "#6b7280"


def apply_figure_style(sizes=(BASE, SECOND, TICK)):
    """Set the rcParams ladder. Call once, before creating the figure."""
    base, secondary, tick = sizes
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": base,
        "axes.labelsize": base, "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": 200, "savefig.dpi": 300,
        "axes.titlelocation": "left",
        "lines.linewidth": 1.2, "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42})


def set_frame(ax):
    for s, v in zip(("top", "right", "bottom", "left"),
                    (False, False, True, True)):
        ax.spines[s].set_visible(v)
    ax.tick_params(direction="out", length=2.5, width=0.6)


def panel_letter(ax, letter, dx=-0.16, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontweight="bold",
            fontsize=LETTER, va="bottom", ha="left")


def effective_pt(fig, authored_pt, doc_width_in=DOC_W):
    """Point size this label will actually render at in the manuscript."""
    return authored_pt * doc_width_in / fig.get_size_inches()[0]


def check_legibility(fig, floor_pt=5.0, doc_width_in=DOC_W):
    """Return every text object that lands below `floor_pt` on the page.

    A figure passes when this returns an empty list. The floor is deliberately
    conservative: 5 pt is roughly the smallest type that survives print
    reduction and PDF rasterization in a two-column layout.
    """
    scale = doc_width_in / fig.get_size_inches()[0]
    bad = []
    for t in fig.findobj(mpl.text.Text):
        if not (t.get_text().strip() and t.get_visible()):
            continue
        eff = t.get_fontsize() * scale
        if eff < floor_pt:
            bad.append((t.get_text()[:40], round(t.get_fontsize(), 1),
                        round(eff, 2)))
    return bad


def check_overlaps(fig):
    """figure-style Sec 9.1 geometric check: visible text boxes must not collide."""
    # Axis label and tick positions are finalized during draw, in the display
    # coordinates of whatever renderer drew them. A preceding savefig() draws
    # at savefig.dpi (300); querying afterwards against the canvas renderer
    # (100 dpi) returns positions from the wrong coordinate system and
    # manufactures overlaps that are not in the image. Redraw first so the
    # positions and the query share one renderer.
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    # Legend entries and tick labels are laid out by matplotlib itself; their
    # window extents include padding that reads as overlap without being one.
    exempt = set()
    ticks = set()
    for ax in fig.axes:
        ticks |= set(ax.get_xticklabels(which="both"))
        ticks |= set(ax.get_yticklabels(which="both"))
        lg = ax.get_legend()
        if lg is not None:
            exempt |= set(lg.findobj(mpl.text.Text))
    for lg in getattr(fig, "legends", []):
        exempt |= set(lg.findobj(mpl.text.Text))

    def shrink(bb, pad=1.0):
        # Text extents carry ~1px of side bearing that touches without the
        # glyphs colliding. Use a small ABSOLUTE pad: a proportional one
        # scales with the label and silently hides real overlaps between
        # long labels (a 10% pad on two 50px labels forgives 10px of
        # collision).
        return mpl.transforms.Bbox.from_extents(bb.x0 + pad, bb.y0,
                                                bb.x1 - pad, bb.y1)

    def extent(t):
        # Annotation.get_window_extent() returns the union of the label and
        # its leader line, so an annotation with a long arrow reports a box
        # many times the size of its text and collides with everything the
        # arrow passes over. Measure the glyphs only.
        if isinstance(t, mpl.text.Annotation):
            return mpl.text.Text.get_window_extent(t, r)
        return t.get_window_extent(r)

    texts = [(t, extent(t)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible() and t not in exempt]
    out = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            # Tick labels vs non-tick text: matplotlib reserves the tick strip,
            # so only flag tick-against-tick collisions (genuinely crowded
            # category axes) and non-tick against non-tick.
            if (a in ticks) != (b in ticks):
                continue
            if shrink(ba).overlaps(shrink(bb)):
                out.append((a.get_text()[:30], b.get_text()[:30]))
    return out
