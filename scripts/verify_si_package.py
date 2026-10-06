#!/usr/bin/env python3
"""Verify the assembled supporting-information package against its sources.

WHY THIS EXISTS
---------------
The submission SI directory is assembled by hand. Nothing regenerated it and
nothing checked it, and in v8.0 that cost four of the five supporting figures:
the manuscript renumbered S3-S5 in a cyclic shift (call-relationships moved
from S5 to S3, comparators from S3 to S4, transformation from S4 to S5), the
figure generator was updated to match, and the .tif files in the package were
never re-copied. The result was a package whose S5_Fig.tif was byte-identical
to the current S3_Fig.tif -- each figure sitting under the wrong number, with
its own README describing something else. Every text gate passed throughout,
because none of them reads a .tif.

This checks the three couplings that failure broke:

  1. every S*_Fig.tif in the package is byte-identical to the render the
     generator currently writes;
  2. the caption title recorded for each SI item in the package README, in
     SI_INDEX.csv and in the figure-mapping record matches the manuscript's
     own @si caption block, which is the authority for numbering;
  3. the file sizes the README advertises match the files actually shipped
     (a cheap tell that a file was replaced without updating its entry).

Exit non-zero on any mismatch. Run with --list to print what it resolved.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
SRC = REPO / "docs" / "manuscript" / "COMAVI_manuscript_source.txt"
GEN = REPO / "figures" / "submission"
MAPPING = REPO / "reference_outputs" / "figure_mapping" / "COMAVI_figure_mapping.csv"

# The package directory is versioned in its name, so resolve it by glob rather
# than hardcoding: a stale literal here would reintroduce exactly the class of
# silent no-op this script exists to prevent (a loop that matches nothing and
# reports success).
def si_dir() -> pathlib.Path:
    hits = sorted((REPO / "submission").glob("COMAVI_*_SI"))
    if len(hits) != 1:
        raise SystemExit(
            "expected exactly one submission/COMAVI_*_SI directory, found %d: %s"
            % (len(hits), [p.name for p in hits]))
    return hits[0]


def sha(p: pathlib.Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manuscript_captions() -> dict[str, str]:
    """{'S1 Fig': 'Stress tests ...'} from the manuscript's @si block."""
    text = SRC.read_text()
    out = {}
    for item, title in re.findall(
            r"@si\s+(S\d+\s+(?:Fig|Table|Text))\s*\|\s*([^|]+)\|", text):
        out[item.strip()] = re.sub(r"\s+", " ", title).strip()
    if not out:
        raise SystemExit("no @si captions found in %s" % SRC)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    SI = si_dir()
    caps = manuscript_captions()
    figs = sorted(k for k in caps if k.endswith("Fig"))
    fail: list[str] = []

    # 1. shipped renders are the current renders
    for item in figs:
        stem = item.replace(" ", "_")
        g, s = GEN / ("%s.tif" % stem), SI / ("%s.tif" % stem)
        if not g.exists():
            fail.append("%s: generator has not written %s" % (item, g.name))
            continue
        if not s.exists():
            fail.append("%s: missing from the package" % item)
            continue
        if sha(g) != sha(s):
            fail.append("%s: package copy differs from the current render "
                        "(regenerate, then copy %s -> %s)" % (item, g, s))
        if args.list:
            print("  %-8s %s  %s" % (item, sha(g)[:16], caps[item][:58]))

    # 2. recorded caption titles agree with the manuscript
    mapped = {}
    with MAPPING.open() as fh:
        for row in csv.DictReader(fh):
            mapped[row["manuscript_item"].strip()] = row["caption_title"].strip()
    for item, title in caps.items():
        if item in mapped and mapped[item] != title:
            fail.append("figure mapping for %s reads %r, manuscript says %r"
                        % (item, mapped[item][:60], title[:60]))

    index = SI / "SI_INDEX.csv"
    if index.exists():
        with index.open() as fh:
            for row in csv.DictReader(fh):
                item = row.get("item", "").strip()
                got = (row.get("caption_title") or "").strip()
                if item in caps and got and not caps[item].startswith(got[:40]):
                    fail.append("SI_INDEX caption for %s reads %r, manuscript "
                                "says %r" % (item, got[:60], caps[item][:60]))

    readme = SI / "README.txt"
    if readme.exists():
        rt = readme.read_text()
        for item in figs:
            stem = item.replace(" ", "_")
            m = re.search(r"%s\.tif[^\n]*?\((\d+) KB\)" % stem, rt)
            if not m:
                fail.append("README has no size entry for %s.tif" % stem)
                continue
            actual = round((SI / ("%s.tif" % stem)).stat().st_size / 1024)
            if int(m.group(1)) != actual:
                fail.append("README says %s.tif is %s KB; it is %d KB"
                            % (stem, m.group(1), actual))
            # the manuscript title must appear beside the filename
            line = next((l for l in rt.splitlines() if "%s.tif" % stem in l), "")
            if caps[item][:38] not in line:
                fail.append("README entry for %s does not carry its "
                            "manuscript caption title" % item)

    if fail:
        print("FAIL: supporting-information package disagrees with its sources")
        for f in fail:
            print("   " + f)
        return 1
    print("PASS: SI package -- %d figures byte-identical to their current "
          "renders; caption titles agree across the manuscript, the figure "
          "mapping, SI_INDEX and README; advertised sizes match"
          % len(figs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
