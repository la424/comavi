#!/usr/bin/env python3
"""Assert CITATION.cff names the newest release this repository documents.

WHY THIS EXISTS
---------------
CITATION.cff carried version "3.0-methods" while v3.0.1-methods was the published
release, and nothing read the file, so nothing could notice. That is the same
shape as every other defect this repository gates against: a published value with
no generator and no check.

The newest release is taken from the RELEASE_NOTES_*.md set rather than from the
git tags, because the notes are what ship inside the archived tarball and are
themselves gated by verify_release_notes_values.py. A tag can be created without
notes; notes without a tag are a release in preparation, and this gate is
satisfied by either as long as the two agree with CITATION.cff.

Known limitation, stated rather than hidden: a version-stamped file inside a
release can only name that release if it is committed BEFORE the tag. The
v3.0.1-methods tarball therefore still contains "3.0-methods"; the archive's own
version metadata, which Zenodo takes from the release, is correct. This gate
makes the lag impossible to repeat silently from the next release onward.

Usage:
    python scripts/verify_citation_version.py
"""
import glob
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]


def _key(v):
    """Sort '3.0.1-methods' after '3.0-methods' and '2.1-methods'."""
    nums = [int(x) for x in re.findall(r"\d+", v.split("-")[0])]
    return nums + [0] * (4 - len(nums))


def main() -> int:
    notes = [pathlib.Path(p).name for p in glob.glob(str(REPO / "RELEASE_NOTES_*.md"))]
    versions = [re.sub(r"^RELEASE_NOTES_v?|\.md$", "", n) for n in notes]
    if not versions:
        print("FAIL: no RELEASE_NOTES_*.md found, so there is nothing to check "
              "CITATION.cff against; a zero-file glob must not report success")
        return 1
    newest = sorted(versions, key=_key)[-1]

    cff = (REPO / "CITATION.cff").read_text()
    m = re.search(r'(?m)^version:\s*"?([^"\n]+)"?\s*$', cff)
    if not m:
        print("FAIL: CITATION.cff has no version field")
        return 1
    declared = m.group(1).strip()

    if declared != newest:
        print("FAIL: CITATION.cff declares version %r but the newest release "
              "notes are for %r" % (declared, newest))
        print("      notes present: %s" % ", ".join(sorted(versions, key=_key)))
        print("      Update CITATION.cff before tagging, or add the notes for the "
              "release it names.")
        return 1
    print("PASS: CITATION.cff version %r matches the newest release notes "
          "(%d release(s) documented)" % (declared, len(versions)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
