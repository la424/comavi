#!/usr/bin/env python3
"""Render the paper of record (.docx) to plain text for claim verification.

The manuscript is hand-authored in Word and is the authoritative source for every
numeric claim. Verifiers need greppable text, so this emits a deterministic
rendering committed alongside the .docx. Regenerate whenever the .docx changes.

Paragraphs and table cells are both emitted, because a substantial number of
manuscript claims live in table captions and cells rather than body prose.
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import sys

from docx import Document

REPO = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_DOCX = REPO / "docs" / "COMAVI_manuscript_current.docx"
DEFAULT_TXT = REPO / "docs" / "COMAVI_manuscript_current.txt"


def render(docx_path: pathlib.Path) -> str:
    doc = Document(str(docx_path))
    out: list[str] = []
    for para in doc.paragraphs:
        text = para.text.strip()
        if text:
            out.append(text)
    for idx, table in enumerate(doc.tables, start=1):
        out.append(f"[TABLE {idx}]")
        for row in table.rows:
            cells = [c.text.strip().replace("\n", " ") for c in row.cells]
            out.append("\t".join(cells))
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--docx", type=pathlib.Path, default=DEFAULT_DOCX)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_TXT)
    ap.add_argument("--check", action="store_true",
                    help="Fail if the rendering differs from the committed text "
                         "instead of rewriting it.")
    args = ap.parse_args()

    if not args.docx.is_file():
        print(f"[FAIL] manuscript not found: {args.docx}", file=sys.stderr)
        return 1

    text = render(args.docx)
    # The supplementary notes and tables moved OUT of the manuscript into separate
    # SI files for PLOS submission, but the claims in them are still published
    # claims and every gate that reads this rendering must still see them.
    # Appending the SI appendix documents keeps the audit's contract intact:
    # every number in the SUBMITTED material traces to a generator, not just
    # every number that happens to sit in the manuscript file.
    si_dir = REPO / "submission" / "COMAVI_v7.12_SI"
    for si in sorted(si_dir.glob("S*_Text.docx")):
        text += "\n[SI %s]\n" % si.stem + render(si)
    for si in sorted(si_dir.glob("S*_Table.xlsx")):
        try:
            import openpyxl
        except ImportError:
            break
        wb = openpyxl.load_workbook(si, read_only=True, data_only=True)
        text += "\n[SI %s]\n" % si.stem
        for ws in wb.worksheets:
            text += "[SHEET %s]\n" % ws.title
            for row in ws.iter_rows(values_only=True):
                text += "\t".join("" if c is None else str(c) for c in row) + "\n"
        wb.close()
    digest = hashlib.sha256(text.encode()).hexdigest()

    if args.check:
        if not args.out.is_file():
            print(f"[FAIL] {args.out} missing; run without --check", file=sys.stderr)
            return 1
        current = args.out.read_text()
        if current != text:
            print(f"[FAIL] {args.out} is stale relative to {args.docx.name}; "
                  f"re-run scripts/extract_manuscript_text.py", file=sys.stderr)
            return 1
        print(f"[OK] {args.out.name} matches {args.docx.name} (sha256 {digest[:16]})")
        return 0

    args.out.write_text(text)
    print(f"wrote {args.out.relative_to(REPO)} "
          f"({len(text.splitlines())} lines, sha256 {digest[:16]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
