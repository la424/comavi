# Public release: what was withheld, and how to recover it

## State
`main` is 1 commit ahead of `origin/main` (e4e0f9b), at `112ca0c`, fast-forward,
working tree clean.
Verified at the current tip: 69 commits scanned, **no withheld path in any
reachable tree**.

The 31st commit (`112ca0c`) catches a defect class no numeric audit can see: a
commit hash cited in a deliverable that a reader cannot fetch. The deck footer
claimed "every number regenerable from the frozen pipeline" next to a commit a
history rewrite had orphaned 26 commits earlier -- unreachable from any ref, and
`git fetch` on it fails. Because no *number* changes when that happens, the
entire suite stayed green. The deck builder is gitignored, so the dead hash was
never public; it reached the built `.pptx`, which is what gets presented.

The footer is now derived by `git rev-parse` at build time and labels its own
publication state: `-dirty` for an unclean tree, `(unpublished)` when HEAD is
not yet an ancestor of `origin/main`. Until this commit is pushed the deck reads
`repo 112ca0c (unpublished)`; rebuild after pushing to clear the marker.

`scripts/verify_provenance_pins.py` gates the class going forward. Four controls
are caught: the dead hash reintroduced, a different orphan pinned, a
one-character typo in the Zenodo tag pin, and a published pin swapped for a
local-only one. Two fixes came from those controls rather than from review -- a
mistyped hash was silently skipped by the first version, and `actions/checkout`
defaults to `depth=1`, which would have made the check skip in CI while passing
locally (`verify.yml` now sets `fetch-depth: 0`).

Separately audited: 20 orphaned commits in the local object store carry withheld
manuscript paths. None was ever a pushed tip, and no pushed tip's tree ever
contained a withheld path, across all 30 tips in the `origin/main` reflog.

The 30th commit (`e4e0f9b`) closes a coverage hole rather than fixing a number.
The README restates every headline the paper reports and no audit read it, so a
value corrected in the manuscript could go stale in the file most readers see
first. `scripts/audit_readme_claims.py` now gates 15 README literals, every one
read from `reference_outputs/COMAVI_numbers_ledger.json` rather than typed, and
matched through the boundary-aware matcher so an over-precise value fails. It
also forbids the superseded 44-set gradient and OR = 6.48.

All 15 gates were mutation-tested and all 15 bind; six substantive controls fail
correctly; crippling the matcher exposes 9 non-binding gates. `verify_audit_binding`
gains a README arm (44 -> 59 gates) that runs on a PUBLIC clone, which previously
got only a skip. `requirements.txt` declares `statsmodels`, needed to regenerate
the BH-adjusted p-values the manuscript tabulates. No reported number changed.

An external review called the README's 99/130 and 0.72 stale "retired" metrics.
They are not: both denominators reconcile to the ledger and to the manuscript's
stated convention. Nothing was retracted.

The 29th commit (`5428f0d`) is the largest correctness fix in this series, and
it is a fix to the *checks*, not the results. No reported number changed.

Both manuscript audits tested prose with substring containment, which cannot see
a value stated to more precision than any analysis produced -- `"0.903"` is
inside `"0.9039"`. Mutation-testing every gate found **45 of 45 digit-ending
checks non-binding**, including both discrimination AUCs, the exact permutation
p, and every screen confidence bound. All 45 passed on corrupted prose.

This was the third defect of one species here, so the fix is structural: one
shared matcher (`scripts/audit_match.py`) with boundary-aware matching, a
partial-update detector for values stated twice where only one mention was
corrected, and a predicate that refuses needles too short to constrain anything.
That predicate immediately rejected three gates that matched almost any text.

It also closed a coverage hole: no audit read the Abstract, so headline values
restated there were gated nowhere. Reversing the Abstract's zero-cell claim with
every number intact now fails.

`scripts/verify_audit_binding.py` makes the sweep permanent in CI, recovering
each audit's literals from the audit itself so the test cannot drift from what it
tests. 44 gates, 0 non-binding; disabling the guards reports 48 and exits 1.

The 28th commit (`edc4210`) makes the bootstrap identity self-testing. The prior
commit claimed algebraic identity in a source comment; a comment is not a test,
so 40 resamples are now recomputed the slow way and must match. The first version
of that check was itself worthless — it re-derived the formula in its own
expression, so breaking the production line left it passing — and a negative
control caught it. Both paths now route through one function. Four controls fail
correctly, including a swapped column pair and a wrong `reindex` fill, which
produce plausible numbers rather than crashes. Output bit-identical.

The 27th commit (`7666a14`) unifies two duplicated estimators. Two statistics
were each computed by two scripts and shipped with two different values, and the
manuscript quoted the wrong one both times because neither literal was gated
against prose. The within-system permutation p now comes only from the exact
enumeration (0.0095, not the Monte-Carlo 0.0097). The specificity-gain bootstrap
was rewritten in closed form — the gain is `(f - fs)/n` over resampled silent
rows — which made a seed-spread sweep cheap and showed the old 4,000-draw upper
bound carried 0.017 of noise; draws raised to 200,000 and reportable precision
cut to `[0.00, +0.24]` with 11% no-gain. Both paired AUCs (0.903 / 0.769) are now
gated in the manuscript audit, not just the deck's.

The 26th commit (`db51cc1`) responds to an external review of the release plan:
CI on the re-deriving track, ablation demotions split by ground-truth class, a
parsed-count audit gate on the demotion sentence, and `scipy` added to
`requirements.txt` (7 scripts import it, including two verifiers — a clean clone
following the README would have failed on ImportError).

## Push
The sandbox cannot reach the macOS keychain, so run this yourself:

    cd ~/mavis_release && git push origin main

Push the branch only. Do **not** use `--all`, `--tags`, or `--mirror`.

## Withheld from public history (18 paths)
Manuscripts `COMAVI_manuscript_v20.md`, `docs/COMAVI_manuscript_v21.{md,docx}`,
`docs/COMAVI_manuscript_v22.{md,docx}`; `docs/COMAVI_outline_v21.md`; the ten
`docs/COMAVI_section_*.md` fragments (including
`COMAVI_section_threshold_operating_points.md`, which was added *and deleted*
inside the unpushed range and would still have been published);
`docs/COMAVI_lab_meeting.pptx`; `scripts/build_lab_meeting_deck.py` (carries 108
manuscript passages as string literals); `scripts/audit_deck_numbers.py`.

The section files were 20–83% verbatim v22 prose, so stripping only the
version-numbered manuscripts would have published most of the paper in pieces.

## Still public, deliberately
Pipeline, canonical outputs, figures, ledger, design decisions, results synthesis,
results specification, writing guide, denominator reconciliation (methods
documentation a reader needs for the 130-vs-131 distinction), and
`scripts/build_manuscript_docx.py` (a general markdown→docx converter).

Residual manuscript prose in public files: one sentence, in
`scripts/analyze_evidence_stratified_agreement.py`'s docstring.

## Clean-clone behaviour (tested)
- `verify_stage6.py` (documented invocation): **11/11**
- `verify_stage6.py --canonical` (re-deriving track): **41/41**
- `verify_denominators`, `verify_tier_construction`, `verify_v30_reconciliation`: PASS
- `audit_evidence_claims`, `audit_tier_energy_claims`: `[SKIP]`, exit 0 — they read
  manuscript prose, so they explain their absence rather than failing.
- CI drift step (`git diff --quiet -- reference_outputs/`): no drift.

`.github/workflows/verify.yml` runs exactly this sequence on every push and PR.
The `--canonical` track is mandatory there, not optional: the default 11-check
track compares stored columns against hardcoded expectations and so passes even
when the shipped scorer can no longer regenerate the stored table — the exact
failure mode of the tier-gating bug. The workflow needs no FoldX.

## Recovery
Files: `_work/manuscript_private/` (17 verbatim copies).
Full pre-rewrite history: `_work/pre-strip-history.bundle` (verified complete,
tip 98a71af, 14 withheld files present).

    git clone _work/pre-strip-history.bundle /tmp/recovered
    # or, into this repo:
    git fetch _work/pre-strip-history.bundle 'refs/heads/main:refs/heads/pre-strip'

Both live under `_work/`, which is gitignored — they cannot be pushed.
The filter-branch backup ref and the `pre-strip-backup` tag were deliberately
deleted: they held the manuscript blobs and a habitual `git push --tags` would
have published them. The bundle is the recovery path.

## Note
`.gitignore` now covers all 18 paths, so a future `git add -A` cannot re-add them.

## Commit 2 of this stack: `add22b6` — restricted-view quantifier disambiguation

`verify_tier_construction.py` emitted the evidence-restricted tier-screen views
under keys (`E1_E3_only`, `direct_or_coupled_only`) whose names read as "every
committed axis qualifies" while the code computed "at least one axis
qualifies". The two differ by 15 variants:

| view | ANY axis | ALL axes |
|---|---|---|
| E1–E3 evidence | n = 42, spec 0.560 | n = 27, spec 0.600 |
| direct or coupled | n = 40, spec 0.565 | n = 25, spec 0.615 |

Both are now emitted under `any_axis_*` / `all_axes_*` with an explicit
`quantifier` field; the old keys remain as ANY aliases. `audit_tier_energy_claims.py`
gates whichever reading the audited section's prose commits to.

A third view seen in draft prose — "fully context-representable variants,
n = 45, specificity 0.571" — is **not computable** from shipped data: no
context-representability column exists in `COMAVI_evidence_ledger.csv` or
`scored_61var_canonical.csv`. Those literals are now forbidden.

### Why this matters for the circulated draft

The draft in circulation states the ALL-quantifier numbers and states them
**correctly** (n = 27/0.600, n = 25/0.615 — both reproduce exactly). It also
states the ungatable context-representable view. Nothing was wrong with the
science; the repo simply could not check it.
