#!/usr/bin/env python3
"""Replace the figures embedded in the manuscript with current generator output.

WHY THIS EXISTS
---------------
Regenerating a figure writes a PNG to figures/; it does not touch the copy
stored inside the .docx. Those two drifted apart silently through an entire
review because no gate read the embedded images (see
scripts/verify_embedded_figures.py, which now does).

WHAT IT DOES
------------
Pairs the inline drawings with the figure captions in document order, swaps in
the bytes of each figure's current file, and rescales the DISPLAY box so the
page layout is preserved: the displayed width is kept exactly and the height is
recomputed from the new image's aspect ratio. Both the wp:extent on the inline
and the a:ext inside the picture's transform are updated -- writers disagree
about which one governs, and leaving them inconsistent distorts the image in
some viewers.

Everything else in the package is copied through byte-for-byte, so styles,
comments, numbering and revision marks survive.

Usage
-----
    PYTHONPATH=. python scripts/embed_figures.py                    # in place
    PYTHONPATH=. python scripts/embed_figures.py --out new.docx
    PYTHONPATH=. python scripts/embed_figures.py --dry-run
"""
import argparse
import io
import pathlib
import re
import shutil
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
DOCX = REPO / "docs" / "COMAVI_manuscript_current.docx"
MAPPING = REPO / "reference_outputs" / "figure_mapping" / "COMAVI_figure_mapping.csv"

NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
}
for k, v in NS.items():
    ET.register_namespace("" if k == "rel" else k, v)
CAPTION_RE = re.compile(r"^\s*((?:Fig \d+|S\d+ Fig))\.\s+\S")


def caption_labels(docx):
    from docx import Document
    d = Document(str(docx))
    out = []
    for el in list(d.element.body):
        if not el.tag.endswith("}p"):
            continue
        txt = "".join(n.text or "" for n in el.iter() if n.tag.endswith("}t"))
        m = CAPTION_RE.match(txt)
        if m and len(txt) > 40:
            out.append(m.group(1))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--docx", type=pathlib.Path, default=DOCX)
    ap.add_argument("--mapping", type=pathlib.Path, default=MAPPING)
    ap.add_argument("--out", type=pathlib.Path, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    from PIL import Image

    mapping = pd.read_csv(args.mapping)
    labels = caption_labels(args.docx)
    zin = zipfile.ZipFile(args.docx)

    rels = ET.fromstring(zin.read("word/_rels/document.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels.iter("{%s}Relationship" % NS["rel"])
              if "image" in (r.get("Type") or "")}

    root = ET.fromstring(zin.read("word/document.xml"))
    inlines = list(root.iter("{%s}inline" % NS["wp"]))
    assert len(inlines) == len(labels) == len(mapping), (
        "%d inline drawings, %d captions, %d mapping rows -- cannot pair"
        % (len(inlines), len(labels), len(mapping)))
    assert labels == list(mapping.manuscript_item), (
        "caption order %s does not match mapping order %s" % (labels, list(mapping.manuscript_item)))

    replace, plan = {}, []
    for inline, label in zip(inlines, labels):
        row = mapping[mapping.manuscript_item == label].iloc[0]
        src = REPO / row.file
        assert src.exists(), "%s: %s is missing" % (label, row.file)

        blip = inline.find(".//{%s}blip" % NS["a"])
        rid = blip.get("{%s}embed" % NS["r"])
        part = "word/" + target[rid].lstrip("/")

        new_bytes = src.read_bytes()
        nw, nh = Image.open(io.BytesIO(new_bytes)).size
        ext = inline.find("{%s}extent" % NS["wp"])
        cx = int(ext.get("cx"))
        new_cy = int(round(cx * nh / nw))
        old_cy = int(ext.get("cy"))
        plan.append((label, pathlib.Path(row.file).name, old_cy, new_cy))

        if not args.dry_run:
            ext.set("cy", str(new_cy))
            aext = inline.find(".//{%s}xfrm/{%s}ext" % (NS["a"], NS["a"]))
            if aext is not None:
                aext.set("cx", str(cx))
                aext.set("cy", str(new_cy))
            replace[part] = new_bytes

    for label, fn, old_cy, new_cy in plan:
        delta = "" if old_cy == new_cy else "  height %d -> %d EMU (%+.1f%%)" % (
            old_cy, new_cy, 100 * (new_cy - old_cy) / old_cy)
        print("  %-8s <- %-46s%s" % (label, fn, delta))

    if args.dry_run:
        print("\ndry run: %d figures would be re-embedded" % len(plan))
        return 0

    out = args.out or args.docx
    tmp = out.with_suffix(".tmp.docx")
    doc_xml = ET.tostring(root, xml_declaration=True, encoding="UTF-8")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename == "word/document.xml":
                zout.writestr(item, doc_xml)
            elif item.filename in replace:
                zout.writestr(item, replace[item.filename])
            else:
                zout.writestr(item, zin.read(item.filename))
    zin.close()
    shutil.move(str(tmp), str(out))
    print("\nwrote %s (%d figures re-embedded)" % (out.relative_to(REPO) if out.is_relative_to(REPO) else out, len(replace)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
