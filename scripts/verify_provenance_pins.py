"""Verify every commit hash cited in a deliverable is fetchable by a reader.

A manuscript that pins its code to a commit hash makes a reproducibility
promise: clone the repo, check out that hash, rerun.  The promise breaks
silently.  Nothing in git warns you that a hash you wrote into prose has since
been dropped by a history rewrite, and no number changes, so every numeric
audit still passes.

This is not hypothetical.  The lab-meeting deck footer cited ``aec59d3`` --
"every number regenerable from the frozen pipeline" -- for 26 commits after a
rewrite had removed that commit from history.  It was unreachable from any
ref.  A reader following the one string offered for reproduction would have got
``fatal: couldn't find remote ref``.

What counts as fetchable is deliberately strict: reachable from
``origin/main``.  Reachable-from-HEAD is not enough, because a local commit is
invisible to a reader until it is pushed.  Tags are accepted when they are
themselves reachable from the published branch.

Scope: the private manuscript sources plus every shipped script.  Figure and
artifact UUIDs are excluded -- they are 32-hex with dashes and are not commits.

Usage:  python scripts/verify_provenance_pins.py
Exit 0 when every cited hash resolves and is published.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Files that make provenance claims to a reader.
TARGETS = [
    "README.md",
    "docs/COMAVI_manuscript_v22.md",
    "scripts/build_lab_meeting_deck.py",
    "scripts/build_manuscript_docx.py",
]

# A hash-shaped token: 7-40 hex chars.  Two exclusions matter:
#   - UUID fragments, which appear in {{artifact:...}} figure markers and are
#     hex but delimited by '-'; requiring a non-'-' neighbour drops them.
#   - decimal numbers, dropped by requiring at least one non-digit hex letter.
HASH = re.compile(r"(?<![0-9a-f-])([0-9a-f]{7,40})(?![0-9a-f-])")

# Only hashes in a PROVENANCE context are citations.  This scoping is what
# lets an unresolvable hash be a failure rather than a silent skip: the repo
# also contains md5/sha1 digests, which are hash-shaped, are not commits, and
# must not be flagged.  Without the context test the only safe behaviour would
# be to ignore every unresolvable token -- which is exactly how a one-character
# typo in a cited hash would go unnoticed.
PROV = re.compile(r"(?i)\b(commit|repo|repository|tag|release|checkout|rev)\b")
CONTEXT = 80  # chars before the hash searched for a provenance keyword


def git(*args: str) -> tuple[int, str]:
    r = subprocess.run(("git", *args), cwd=REPO, capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def is_commit(h: str) -> bool:
    rc, out = git("cat-file", "-t", h)
    return rc == 0 and out == "commit"


def published(h: str) -> bool:
    """True when a reader who clones the public repo gets this commit."""
    return git("merge-base", "--is-ancestor", h, "origin/main")[0] == 0


def cited_hashes(text: str) -> set[str]:
    """Hashes presented to a reader as a commit citation."""
    out = set()
    for m in HASH.finditer(text):
        h = m.group(1)
        if h.isdigit():          # pure decimal, e.g. a DOI fragment
            continue
        if not PROV.search(text[max(0, m.start() - CONTEXT):m.start()]):
            continue             # a digest, not a provenance claim
        out.add(h)
    return out


def main() -> int:
    rc, _ = git("rev-parse", "--git-dir")
    if rc != 0:
        print("verify_provenance_pins: not a git checkout, skipping")
        return 0
    if git("rev-parse", "--verify", "origin/main")[0] != 0:
        print("verify_provenance_pins: no origin/main ref, skipping")
        return 0
    # A shallow clone (actions/checkout defaults to depth=1) has no ancestry,
    # so every cited hash would look unresolvable and this would fail for a
    # reason that has nothing to do with the citations.  Skip loudly: a silent
    # pass here would be worse than the defect the script exists to catch.
    if git("rev-parse", "--is-shallow-repository")[1] == "true":
        print("verify_provenance_pins: SHALLOW clone -- ancestry unavailable, "
              "skipping.  CI must set fetch-depth: 0 for this check to run.",
              file=sys.stderr)
        return 0

    failures: list[str] = []
    checked = 0
    for rel in TARGETS:
        p = REPO / rel
        if not p.exists():       # manuscript is withheld from public clones
            continue
        for h in sorted(cited_hashes(p.read_text())):
            checked += 1
            if not is_commit(h):
                failures.append(
                    f"{rel}: cites {h}, which resolves to no commit object "
                    f"in this repo -- a mistyped or rewritten hash"
                )
            elif not published(h):
                failures.append(
                    f"{rel}: cites {h}, which is NOT reachable from "
                    f"origin/main -- a reader cannot fetch it"
                )

    print(f"provenance pins : {checked} commit citations checked")
    for f in failures:
        print(f"  FAIL  {f}", file=sys.stderr)
    if failures:
        print(f"\n{len(failures)} unfetchable commit citation(s)", file=sys.stderr)
        return 1
    print("all cited commits are reachable from origin/main")
    return 0


if __name__ == "__main__":
    sys.exit(main())
