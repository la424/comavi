#!/usr/bin/env python3
"""Do the figures INSIDE the manuscript match the figures the generators write?

WHY THIS EXISTS
---------------
Every other gate in this repository reads the manuscript's TEXT -- paragraphs
and table cells. A .docx stores its figures as embedded image parts under
word/media/, which no text audit can see. So the document went through a full
numerical review with corrected prose sitting next to pre-correction images:
every one of the eleven embedded figures was the version from before the v7.6
ledger correction and the v7.8 figure repairs. Two of them were caught by a
human reading the document, not by any check here.

Regenerating a PNG on disk does not update the copy inside the document. This
gate is the only thing that notices.

METHOD
------
Images are read in document order (following the drawing relationships, not the
media/ filenames, which are assigned by the writer and carry no meaning) and
paired positionally with the figure captions, also in document order. That
pairing is asserted to be one-to-one and the caption labels are asserted to
match the expected manuscript sequence, so a reordered or dropped figure fails
here rather than silently shifting the comparison.

Comparison is on DECODED PIXELS, not file bytes. PNG is lossless, so a writer
that re-containers an image leaves the pixels identical and the gate stays
quiet; an image that is actually a different render fails. A size mismatch is
reported separately from a content mismatch because they have different causes
-- a resized embed versus a stale one.

Usage
-----
    PYTHONPATH=. python scripts/verify_embedded_figures.py
    PYTHONPATH=. python scripts/verify_embedded_figures.py --docx path/to.docx
"""
import argparse
import hashlib
import io
import pathlib
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
DOCX = REPO / "docs" / "COMAVI_manuscript_current.docx"
MAPPING = REPO / "reference_outputs" / "figure_mapping" / "COMAVI_figure_mapping.csv"

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
RELNS = "{http://schemas.openxmlformats.org/package/2006/relationships}"
REMBED = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
CAPTION_RE = re.compile(r"^\s*((?:Fig \d+|S\d+ Fig))\.\s+\S")


def pixel_digest(data):
    """SHA256 over decoded RGB pixels, so lossless re-encoding does not trip.

    Pillow is a hard requirement of this gate, not an optional extra: without a
    decoder there is no way to tell a re-containered image from a different
    render, and a byte comparison would raise a false alarm every time the
    document is opened and saved. Fail loudly rather than silently degrading to
    a weaker check -- a verification script that quietly stops verifying is
    worse than one that stops.
    """
    try:
        from PIL import Image
    except ImportError:
        raise SystemExit(
            "FAIL: this gate needs Pillow to decode the embedded images.\n"
            "  pip install pillow  (it is in requirements.txt)\n"
            "  Refusing to fall back to a byte comparison, which would fail on "
            "any document that has been opened and saved.")
    im = Image.open(io.BytesIO(data)).convert("RGB")
    return hashlib.sha256(im.tobytes()).hexdigest(), im.size


def embedded_images(docx):
    z = zipfile.ZipFile(docx)
    rels = ET.fromstring(z.read("word/_rels/document.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels.iter(f"{RELNS}Relationship")
              if "image" in (r.get("Type") or "")}
    doc = ET.fromstring(z.read("word/document.xml"))
    out = []
    for blip in doc.iter(f"{A}blip"):
        rid = blip.get(f"{REMBED}embed")
        if rid in target:
            part = "word/" + target[rid].lstrip("/")
            out.append((target[rid], z.read(part)))
    return out


def caption_labels(docx):
    from docx import Document
    d = Document(str(docx))
    labels = []
    for el in list(d.element.body):
        if not el.tag.endswith("}p"):
            continue
        txt = "".join(n.text or "" for n in el.iter() if n.tag.endswith("}t"))
        m = CAPTION_RE.match(txt)
        if m and len(txt) > 40:
            labels.append(m.group(1))
    return labels


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--docx", type=pathlib.Path, default=DOCX)
    ap.add_argument("--mapping", type=pathlib.Path, default=MAPPING)
    ap.add_argument("--check", action="store_true",
                    help="accepted for symmetry with the build_* generators")
    args = ap.parse_args()

    if not args.docx.exists():
        print("FAIL: %s not found" % args.docx)
        return 1
    if not args.mapping.exists():
        print("FAIL: figure mapping %s not found" % args.mapping)
        return 1

    mapping = pd.read_csv(args.mapping)
    # PLOS uploads supporting information as separate files and lists only its
    # CAPTIONS at the end of the manuscript, so the S# Fig items are named in the
    # document but are deliberately not embedded in it. Only the main figures are
    # paired against embedded images here; the S# Fig rows are checked for the
    # opposite property, that they are NOT embedded, so a supplementary figure
    # accidentally left inside the manuscript fails rather than passing quietly.
    main_expected = [x for x in mapping.manuscript_item if str(x).startswith("Fig ")]
    si_expected = [x for x in mapping.manuscript_item if not str(x).startswith("Fig ")]
    images = embedded_images(args.docx)
    labels = caption_labels(args.docx)
    main_labels = [x for x in labels if x.startswith("Fig ")]
    si_labels = [x for x in labels if not x.startswith("Fig ")]

    problems = []
    if len(images) != len(main_labels):
        print("FAIL: %d embedded images but %d main-figure captions -- cannot pair them"
              % (len(images), len(main_labels)))
        print("  main-figure captions: %s" % main_labels)
        print("  supplementary captions present (must not be embedded): %s" % si_labels)
        return 1
    if main_labels != main_expected:
        print("FAIL: main-figure caption sequence does not match the mapping")
        print("  captions: %s" % main_labels)
        print("  mapping : %s" % main_expected)
        return 1
    if sorted(si_labels) != sorted(si_expected):
        print("FAIL: supplementary figure captions do not match the mapping")
        print("  captions: %s" % sorted(si_labels))
        print("  mapping : %s" % sorted(si_expected))
        return 1
    for label in si_expected:
        src = REPO / mapping[mapping.manuscript_item == label].iloc[0].file
        if not src.exists():
            problems.append("%s: generator output %s is missing" % (label, src))
    labels = main_labels

    print("pairing %d embedded images with %d captions, in document order\n" % (len(images), len(labels)))
    for (target, data), label in zip(images, labels):
        row = mapping[mapping.manuscript_item == label].iloc[0]
        src = REPO / row.file
        if not src.exists():
            problems.append("%s: generator output %s is missing" % (label, row.file))
            print("  %-8s MISSING  %s" % (label, row.file))
            continue
        emb_hash, emb_size = pixel_digest(data)
        src_hash, src_size = pixel_digest(src.read_bytes())
        if emb_hash == src_hash:
            print("  %-8s OK       %s" % (label, pathlib.Path(row.file).name))
        elif emb_size != src_size:
            problems.append("%s: embedded image is %dx%d, generator writes %dx%d (%s)"
                            % (label, emb_size[0], emb_size[1], src_size[0], src_size[1], row.file))
            print("  %-8s SIZE     embedded %dx%d vs file %dx%d" % (label, *emb_size, *src_size))
        else:
            problems.append("%s: embedded image differs in content from %s -- it is a "
                            "different render, most likely predating a figure correction"
                            % (label, row.file))
            print("  %-8s STALE    same size, different pixels: %s" % (label, pathlib.Path(row.file).name))

    if problems:
        print("\nFAIL: %d embedded figure(s) do not match their generator output" % len(problems))
        for p in problems:
            print("  - %s" % p)
        print("\nRe-embed with scripts/embed_figures.py, then re-run this gate.")
        return 1

    print("\nPASS: all %d embedded figures match the files their generators write" % len(images))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
