#!/usr/bin/env python3
"""Build S1-S3 Text (.docx) from their plain-text sources.

WHY THIS EXISTS
---------------
Until v11 the three supporting texts were hand-edited Word files. Nothing
rebuilt them, so when the records moved they kept the old values: the v9.2 S1
Text still described a 99-expectation, all-row ledger after the benchmark had
become the 56 graded variants, and S3 Text described a correction to KRAS
G12D and G12V that a later, source-verified correction had reversed. The
manuscript itself has been built from a source file since v10; this script
gives the supplement the same footing:

  docs/manuscript/si/S{1,2,3}_Text_source.txt  ->  submission/COMAVI_v*_SI/S{1,2,3}_Text.docx

The sources use the manuscript's markup (scripts/build_manuscript.py), plus a
few directives of their own:

  @title <text>             first line of the file, bold
  @contents                 table of the section headings that follow
  ## <heading>              section heading
  @code <text>              one line in a fixed-width font
  @tablecap <text>          caption above a table ("Table A. ...")
  @tabhead a | b | c        table header row
  @tabrow a | b | c         table body row (consecutive rows form one table)
  @references               the file's own numbered reference list
  [@key; @key]              citation, numbered within the file

Each supporting file carries its own reference list, numbered from 1 in order
of first appearance, because PLOS requires supporting files to be
self-contained. Keys resolve against docs/manuscript/references.json and, for
references cited only in the supplement, docs/manuscript/si/si_references.json.

The values in the sources were written by a script that reads them from the
committed records and asserts the premises the prose relies on; the values and
their record paths are listed in docs/audit/COMAVI_v11_si_text_values.json.

Usage:
    python scripts/build_si_texts.py           # write the three .docx files
    python scripts/build_si_texts.py --check   # rebuild in memory, compare text, write nothing
"""
from __future__ import annotations

import argparse
import io
import json
import pathlib
import re
import sys
import tempfile
import zipfile

from docx import Document
from docx.shared import Pt

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import build_manuscript as bm  # noqa: E402

SI_SRC = REPO / "docs" / "manuscript" / "si"
TEXTS = ("S1_Text", "S2_Text", "S3_Text")
REFS = REPO / "docs" / "manuscript" / "references.json"
SI_REFS = SI_SRC / "si_references.json"
PROVENANCE = SI_SRC / "si_build_provenance.json"
CODE_FONT, CODE_PT = "Courier New", 8


def si_dir() -> pathlib.Path:
    hits = sorted((REPO / "submission").glob("COMAVI_v*_SI"))
    if len(hits) != 1:
        raise SystemExit("expected exactly one submission/COMAVI_v*_SI directory, found %d: %s"
                         % (len(hits), [h.name for h in hits]))
    return hits[0]


def references() -> dict:
    refs = json.loads(REFS.read_text())
    if SI_REFS.exists():
        extra = json.loads(SI_REFS.read_text())
        clash = sorted(set(extra) & set(refs))
        if clash:
            raise KeyError("keys defined in both reference files: %s" % clash)
        refs.update(extra)
    return refs


def number_citations(blocks) -> dict:
    """First-appearance order within one supporting file."""
    order = []
    for kind, payload, _ in blocks:
        if kind in ("p", "@tabrow", "@tablecap", "@tabhead"):
            for m in bm.CITE.finditer(payload):
                for k in bm.cite_keys(m.group(1)):
                    if k not in order:
                        order.append(k)
    return {k: i + 1 for i, k in enumerate(order)}


def cells(payload: str) -> list[str]:
    return [c.strip() for c in payload.split("|")]


def build_one(name: str, out_path: pathlib.Path, refs: dict) -> dict:
    src = SI_SRC / ("%s_source.txt" % name)
    blocks = bm.parse_source(src.read_text())
    numbers = number_citations(blocks)
    unknown = [k for k in numbers if k not in refs]
    if unknown:
        raise KeyError("%s: citation keys missing from the reference files: %s" % (name, unknown))
    if numbers and not any(k == "@references" for k, _, _ in blocks):
        raise ValueError("%s cites references but has no @references block" % name)

    doc = Document(str(bm.TEMPLATE))
    bm.clear_body(doc)
    sections = [p for k, p, _ in blocks if k == "h2"]
    prov, i = [], 0
    while i < len(blocks):
        kind, payload, origin = blocks[i]
        if kind == "@title":
            bm.add_runs(doc.add_paragraph(), payload, bold=True)
        elif kind == "@contents":
            rows = []
            for s in sections:
                label, _, title = s.partition(". ")
                rows.append([label, title])
            bm.write_table(doc, ["Section", "Title"], rows, left_cols=(0, 1))
            doc.add_paragraph()
        elif kind == "h2":
            doc.add_paragraph(payload, style="Heading 2")
        elif kind == "p":
            text = bm.render_citations(payload, numbers)
            bm.add_runs(doc.add_paragraph(), text)
            prov.append({"text": bm.plain(text), "from": origin})
        elif kind == "@code":
            par = doc.add_paragraph()
            par.paragraph_format.space_after = Pt(2)
            run = par.add_run(payload)
            run.font.name, run.font.size = CODE_FONT, Pt(CODE_PT)
        elif kind == "@tablecap":
            bm.add_runs(doc.add_paragraph(), bm.render_citations(payload, numbers), bold=True)
            prov.append({"text": bm.plain(bm.render_citations(payload, numbers)), "from": "caption"})
        elif kind == "@tabhead":
            header = cells(payload)
            rows = []
            while i + 1 < len(blocks) and blocks[i + 1][0] == "@tabrow":
                i += 1
                row = cells(bm.render_citations(blocks[i][1], numbers))
                if len(row) != len(header):
                    raise ValueError("%s: table row has %d cells, header %d: %s"
                                     % (name, len(row), len(header), blocks[i][1][:80]))
                rows.append(row)
            bm.write_table(doc, header, rows, left_cols=tuple(range(len(header))))
            doc.add_paragraph()
        elif kind == "@tabrow":
            raise ValueError("%s: @tabrow without a preceding @tabhead" % name)
        elif kind == "@references":
            doc.add_paragraph("References", style="Heading 2")
            for k, n in sorted(numbers.items(), key=lambda kv: kv[1]):
                bm.add_runs(doc.add_paragraph(), "%d. %s" % (n, refs[k]["text"]))
        else:
            raise ValueError("%s: unknown directive %s" % (name, kind))
        i += 1
    doc.save(str(out_path))
    out_path.write_bytes(normalize_zip(out_path.read_bytes()))
    return {"numbers": numbers, "provenance": prov}


# Fixed zip-entry timestamps and document dates: an unchanged build is then
# byte-identical, so the package checksum manifest (SHA256SUMS.txt, written by
# build_supplement.py) stays valid across rebuilds. Shared with that script.
STAMP = (2026, 1, 1, 0, 0, 0)
STAMP_ISO = "2026-01-01T00:00:00Z"


def normalize_zip(data: bytes) -> bytes:
    src = zipfile.ZipFile(io.BytesIO(data))
    out = io.BytesIO()
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as dst:
        for info in src.infolist():
            blob = src.read(info.filename)
            if info.filename == "docProps/core.xml":
                blob = re.sub(rb"(<dcterms:(created|modified)[^>]*>)[^<]*(</dcterms:\2>)",
                              lambda m: m.group(1) + STAMP_ISO.encode() + m.group(3), blob)
            zi = zipfile.ZipInfo(info.filename, date_time=STAMP)
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.external_attr = 0o644 << 16
            dst.writestr(zi, blob)
    return out.getvalue()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="rebuild and compare text with the committed .docx files; write nothing")
    args = ap.parse_args()
    refs, out_dir = references(), si_dir()
    if args.check:
        bad = []
        with tempfile.TemporaryDirectory() as td:
            for name in TEXTS:
                tmp = pathlib.Path(td) / (name + ".docx")
                build_one(name, tmp, refs)
                committed = out_dir / (name + ".docx")
                if not committed.exists():
                    bad.append("%s missing" % committed.relative_to(REPO))
                elif bm.doc_text(tmp) != bm.doc_text(committed):
                    bad.append("%s differs from a fresh build of %s_source.txt" % (committed.name, name))
        for b in bad:
            print("FAIL: %s" % b)
        if bad:
            return 1
        print("PASS: %s reproduce from their sources" % ", ".join(t + ".docx" for t in TEXTS))
        return 0
    record = {}
    for name in TEXTS:
        record[name] = build_one(name, out_dir / (name + ".docx"), refs)
        print("wrote %s: %d references, %d paragraphs"
              % ((out_dir / (name + ".docx")).relative_to(REPO), len(record[name]["numbers"]),
                 len(record[name]["provenance"])))
    PROVENANCE.write_text(json.dumps(record, ensure_ascii=False, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
