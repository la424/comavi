#!/usr/bin/env python3
"""Do the submission figure files meet the PLOS figure specification?

WHY THIS EXISTS
---------------
Until v7.14 nothing checked a figure file at all. The manuscript's figures
were drawn up to 12.7 in wide and shrunk to the text column, so nearly every
label printed below the journal's 8 pt minimum, in fonts the journal does not
accept, and no gate could see it: every other check here reads text.

WHAT IT CHECKS
--------------
Files (always; no plotting library needed, so CI runs it):
  * every figure the manuscript source cites has a TIFF upload file and a PNG
    embed copy in figures/submission/, and no stray figure file is present
  * TIFF: LZW-compressed, RGB with no alpha channel, 600 dpi, under 10 MB,
    2.63-7.5 in wide and at most 8.75 in tall at that resolution
  * PNG: RGB, 300 dpi, the same printed size as its TIFF to within one pixel
  * figures/submission/SHA256SUMS.txt verifies, so a file cannot change
    without its manifest

Drawing (--render; needs Arial, so it runs where the figures are built):
  * rebuilds each figure in memory and applies the generator's own checks:
    Arial only, every text element 8-12 pt at printed size, no title inside
    the image, no colliding labels

Negative controls (always): synthetic files that break one rule each -- RGBA,
uncompressed, 300 dpi, 8 in wide -- must each fail, and with --render a
7 pt label, a DejaVu Sans label and an axes title must each be flagged. A
check that cannot fail is not a check.

Usage (repo root):
    python scripts/verify_figure_specs.py            # files + negative controls
    python scripts/verify_figure_specs.py --render   # also rebuild and inspect the drawings
"""
import argparse
import hashlib
import io
import pathlib
import re
import sys
import tempfile

from PIL import Image

REPO = pathlib.Path(__file__).resolve().parents[1]
FIGDIR = REPO / "figures" / "submission"
SOURCE = REPO / "docs" / "manuscript" / "COMAVI_manuscript_source.txt"
TIFF_DPI, PNG_DPI = 600, 300
W_MIN, W_MAX, H_MAX = 2.63, 7.5, 8.75
MAX_BYTES = 10 * 1024 * 1024


def cited_figures():
    """File stems for every figure the manuscript source declares."""
    text = SOURCE.read_text()
    main = re.findall(r"^@figure (Fig \d+) \|", text, re.M)
    si = re.findall(r"^@si (S\d+) Fig \|", text, re.M)
    return [m.replace(" ", "") for m in main] + ["%s_Fig" % s for s in si]


def check_tif(path):
    bad = []
    try:
        im = Image.open(path)
    except Exception as exc:  # noqa: BLE001
        return ["%s: not readable (%s)" % (path.name, exc)]
    if im.format != "TIFF":
        bad.append("%s: format %s, not TIFF" % (path.name, im.format))
    if im.mode != "RGB":
        bad.append("%s: mode %s, not RGB without alpha" % (path.name, im.mode))
    comp = im.info.get("compression")
    if comp != "tiff_lzw":
        bad.append("%s: compression %s, not LZW" % (path.name, comp))
    dpi = im.info.get("dpi")
    if not dpi or any(abs(float(d) - TIFF_DPI) > 0.5 for d in dpi):
        bad.append("%s: resolution %s, not %d dpi" % (path.name, dpi, TIFF_DPI))
    w_in, h_in = im.size[0] / TIFF_DPI, im.size[1] / TIFF_DPI
    if not (W_MIN <= w_in <= W_MAX):
        bad.append("%s: %.2f in wide, outside %.2f-%.2f in" % (path.name, w_in, W_MIN, W_MAX))
    if h_in > H_MAX:
        bad.append("%s: %.2f in tall, over %.2f in" % (path.name, h_in, H_MAX))
    if path.stat().st_size >= MAX_BYTES:
        bad.append("%s: %.1f MB, not under 10 MB" % (path.name, path.stat().st_size / 1e6))
    return bad


def check_png(png, tif):
    bad = []
    im = Image.open(png)
    if im.mode != "RGB":
        bad.append("%s: mode %s, not RGB" % (png.name, im.mode))
    dpi = im.info.get("dpi")
    if not dpi or any(abs(float(d) - PNG_DPI) > 0.5 for d in dpi):
        bad.append("%s: resolution %s, not %d dpi" % (png.name, dpi, PNG_DPI))
    if tif.exists():
        t = Image.open(tif)
        for k in (0, 1):
            if abs(im.size[k] / PNG_DPI - t.size[k] / TIFF_DPI) > 1.5 / PNG_DPI:
                bad.append("%s: printed size differs from %s" % (png.name, tif.name))
                break
    return bad


def check_manifest():
    man = FIGDIR / "SHA256SUMS.txt"
    if not man.exists():
        return ["SHA256SUMS.txt missing"]
    bad = []
    listed = set()
    for line in man.read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        listed.add(name.strip())
        p = FIGDIR / name.strip()
        if not p.exists():
            bad.append("manifest lists a missing file: %s" % name)
        elif hashlib.sha256(p.read_bytes()).hexdigest() != digest:
            bad.append("manifest checksum mismatch: %s" % name)
    for p in FIGDIR.glob("*"):
        if p.suffix in (".tif", ".png") and p.name not in listed:
            bad.append("figure file not in the manifest: %s" % p.name)
    return bad


def check_files():
    bad = []
    stems = cited_figures()
    if not stems:
        return ["no figures found in %s" % SOURCE.relative_to(REPO)]
    for stem in stems:
        tif, png = FIGDIR / (stem + ".tif"), FIGDIR / (stem + ".png")
        if not tif.exists():
            bad.append("%s.tif missing" % stem)
        else:
            bad += check_tif(tif)
        if not png.exists():
            bad.append("%s.png missing" % stem)
        else:
            bad += check_png(png, tif)
    stray = sorted(p.name for p in FIGDIR.glob("*") if p.suffix in (".tif", ".png") and p.stem not in stems)
    bad += ["stray figure file not cited by the manuscript: %s" % s for s in stray]
    bad += check_manifest()
    return bad, stems


def negative_controls_files():
    """Each synthetic file breaks one rule and must be caught."""
    failures = []
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        base = Image.new("RGB", (int(6.5 * TIFF_DPI), int(3 * TIFF_DPI)), "white")
        cases = {
            "rgba": (base.convert("RGBA"), dict(compression="tiff_lzw", dpi=(600, 600))),
            "uncompressed": (base, dict(dpi=(600, 600))),
            "300dpi": (base, dict(compression="tiff_lzw", dpi=(300, 300))),
            "too_wide": (Image.new("RGB", (int(8 * TIFF_DPI), 600), "white"),
                         dict(compression="tiff_lzw", dpi=(600, 600))),
        }
        good = td / "good.tif"
        base.save(good, format="TIFF", compression="tiff_lzw", dpi=(600, 600))
        if check_tif(good):
            failures.append("a compliant control file was rejected: %s" % check_tif(good))
        for name, (im, kw) in cases.items():
            p = td / (name + ".tif")
            im.save(p, format="TIFF", **kw)
            if not check_tif(p):
                failures.append("negative control '%s' was NOT caught" % name)
    return failures


def check_render(stems):
    sys.path.insert(0, str(FIGDIR))
    import make_submission_figures as msf  # noqa: E402
    import matplotlib.pyplot as plt  # noqa: E402

    bad = []
    for stem in stems:
        msf.style()
        fig = msf.FIGURES[stem]()
        v = msf.spec_violations(fig)
        plt.close(fig)
        bad += ["%s: %s" % (stem, x) for x in v]
    # Negative controls on the drawing checks.
    failures = []
    for name, build in [
            ("7 pt label", lambda ax: ax.text(0.5, 0.5, "small", fontsize=7)),
            ("DejaVu Sans label", lambda ax: ax.text(0.5, 0.5, "wrong font", family="DejaVu Sans")),
            ("axes title", lambda ax: ax.set_title("A title inside the image")),
            ("colliding labels", lambda ax: (ax.text(0.5, 0.5, "first label"), ax.text(0.51, 0.5, "second"))),
    ]:
        msf.style()
        fig, ax = plt.subplots(figsize=(3.4, 2.5))
        build(ax)
        caught = bool(msf.spec_violations(fig))
        plt.close(fig)
        if not caught:
            failures.append("drawing negative control '%s' was NOT caught" % name)
    return bad, failures


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--render", action="store_true", help="also rebuild and inspect each drawing (needs Arial)")
    args = ap.parse_args()
    bad, stems = check_files()
    controls = negative_controls_files()
    print("figures declared by the manuscript source: %s" % ", ".join(stems))
    if args.render:
        rbad, rcontrols = check_render(stems)
        bad += rbad
        controls += rcontrols
    for x in bad:
        print("  FAIL %s" % x)
    for x in controls:
        print("  CONTROL FAILED %s" % x)
    if bad or controls:
        print("FAIL: %d specification problem(s), %d negative control(s) not caught" % (len(bad), len(controls)))
        return 1
    print("PASS: %d figures meet the file specification%s; all negative controls caught"
          % (len(stems), " and the drawing specification" if args.render else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
