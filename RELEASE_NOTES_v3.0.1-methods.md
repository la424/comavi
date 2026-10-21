# COMAVI v3.0.1-methods

A patch release. **No published value differs from v3.0-methods.** Cite either; the
concept DOI resolves to this one.

It exists because v3.0-methods shipped at a commit that fails one of its own checks.
`verify_no_duplicated_headlines.py` flagged `RELEASE_NOTES_v3.0-methods.md` for
restating the mechanism-pattern agreement, and that gate was right: an unchecked
document that restates a headline is how this project previously published a figures
README quoting a cohort that no longer existed, and a results summary whose template
hardcoded a fraction beside the computed rate it was meant to describe. A release whose
notes claim every value is re-derived in CI should not ship with that check red.

Release notes cannot be forbidden from stating headline numbers, so they are now allowed
by pattern and gated instead:

- `scripts/verify_release_notes_values.py` rebuilds every current-column value from the
  committed records and asserts it appears in the notes, in the bold form the current
  column uses. Needles are computed, never typed. A file is held to the rows it actually
  claims, and every tracked quantity must be claimed by at least one notes file, so a row
  cannot quietly disappear.
- `verify_no_duplicated_headlines.py` gained `ALLOWED_PATTERNS`, so a future release needs
  no edit. Both run in the same CI step: an allowance without its gate is worse than no
  allowance.

Negative controls: a stale current value, a current value left unbolded, a changed prose
interval, and an ungated document restating a headline are each caught.

The previous-release column in the v3.0-methods notes is deliberately not checked against
the records — those values are superseded by construction. They were verified against what
`v2.1-methods` actually published, read out of that tag.

## Everything else in this release

Unchanged from v3.0-methods. See `RELEASE_NOTES_v3.0-methods.md` for the ground-truth
corrections, the column-selection defect that fabricated a published value, the
supplementary evidence stratification, and the verification-infrastructure repairs.
