"""Deterministic figure output. Import before any savefig.

WHY THIS EXISTS
---------------
Matplotlib's PDF backend stamps a CreationDate into every file it writes, so
rebuilding an unchanged figure produces a byte-different PDF whose only
difference is that timestamp. Five of the manuscript's PDF sidecars showed
exactly that: identical byte length, identical content stream, one differing
field, while their PNG siblings were byte-identical.

The cost is not cosmetic. The checksum manifests in figures/ and
figures/isds_v1/ can then never verify after a rebuild, so the repository
loses the property every other generator here maintains -- that re-running it
rewrites no committed output -- for every PDF. Worse, a fresh clone that
regenerates figures sees five mismatches it cannot tell apart from real drift.

Matplotlib honours SOURCE_DATE_EPOCH, the reproducible-builds convention, for
that field. Setting it here fixes every savefig in every generator that
imports this module, instead of editing 23 call sites.

The pinned value is this repository's first commit (2026-06-22), a real date
rather than an invented constant. setdefault, so a build system that already
pins SOURCE_DATE_EPOCH keeps its own value.

USAGE
    import _repro  # noqa: F401  (imported for its side effect)

Generators under figures/src/ can import it directly; those elsewhere need
figures/src on sys.path first, as figures/isds_v1/make_figures.py does.
"""

import os

REPO_FIRST_COMMIT_EPOCH = "1782152076"   # 2026-06-22

os.environ.setdefault("SOURCE_DATE_EPOCH", REPO_FIRST_COMMIT_EPOCH)
