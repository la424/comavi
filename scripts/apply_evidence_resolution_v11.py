#!/usr/bin/env python3
"""Apply the v11 evidence-resolution rule to the evidence ledger and the canonical.

WHY THIS EXISTS
---------------
An expectation can be no more specific than the measurement behind it. A
co-immunoprecipitation, pull-down or two-hybrid assay reports how much complex
forms. A reduced signal can come from weaker binding or from a less stable
assembled complex, and COMAVI computes these two separately (the AnalyseComplex
interaction energy and the BuildModel stability of the assembled complex). When
such an assay is the only evidence for an assembled-complex effect, a
binding-specific expectation asserts more than was measured, and a call in the
measured direction on either assembled axis matches what was measured.

The rule is defined by the evidence, not by any call: it applies to every
expectation that meets it, whether or not COMAVI benefits, and to every graded
strategy, the constant baselines included (apply_concordance_v5 implements the
grading; see EITHER_ASSEMBLED_AXIS there).

It only relaxes over-specific expectations. A normal complex amount rules out a
binding loss but cannot establish that the assembled complex is stable: MLH1
R755W has a PMS2 co-IP at 82% of wild type and a destabilized assembled complex
by MAVE abundance. So no commitment is added on the non-disruptive side; the
four neutral binding expectations that rest on such an assay are recorded in the
new ledger column and left unchanged.

WHAT IT CHANGES
---------------
1. Evidence ledger: a column complex_amount_assay naming the assay for every
   binding expectation whose evidence is a measurement of how much complex
   forms. Each named assay must appear in that row's evidence basis or citation,
   and every binding row whose text names such an assay must be classified here.
2. Canonical: a column assembled_expectation_scope, "either_axis" for a variant
   whose binding expectation is a destabilization or stabilization resting on
   such an assay and whose assembled-complex stability axis was not tested.
3. Grades: every grade is recomputed through the same regrade() as
   scripts/apply_rubric_table3_v11.py, and every change is written to
   docs/audit/COMAVI_v11_evidence_resolution.csv. Any change outside the
   either_axis variants stops the run.
4. Direction agreement: the same expectations are scored once, jointly, on
   either assembled axis ('assembled_either' row of
   apply_concordance_v5.structural_agreement_by_axis), so each axis's own row
   is tested only against expectations specific to that axis. The stored
   per-variant structural_agreement_* and directional_agreement_* columns are
   recomputed at every threshold, and every moved cell is written to
   docs/audit/COMAVI_v11_evidence_resolution_agreement.csv. Any change outside
   the either_axis variants stops the run.

NOT APPLIED
-----------
VHL L158Q commits both assembled axes on "Disrupts VCB complex assembly"
(Stebbins 1999 Science, a crystal-structure paper, closed access). The ledger
names no assay, so the rule is not applied and both commitments stand.

Usage:
    python scripts/apply_evidence_resolution_v11.py           # apply
    python scripts/apply_evidence_resolution_v11.py --check   # verify applied
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import apply_concordance_v5 as ac  # noqa: E402
import apply_rubric_table3_v11 as table3  # noqa: E402

LEDGER = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
AUDIT = REPO / "docs" / "audit" / "COMAVI_v11_evidence_resolution.csv"
# Per-variant agreement cells moved by scoring complex-amount expectations
# jointly in direction agreement (decided after the grading change; see agreement()).
AGREEMENT_AUDIT = REPO / "docs" / "audit" / "COMAVI_v11_evidence_resolution_agreement.csv"
LEDGER_COLUMN = "complex_amount_assay"
SCOPE_COLUMN = "assembled_expectation_scope"

# Binding expectations whose evidence is a measurement of how much complex
# forms. Every "; "-separated name must appear in the row's basis or citation.
COMPLEX_AMOUNT = {
    ("brca1_bard1", "V11G"): "mammalian two-hybrid",
    ("brca1_bard1", "E100Q"): "Y2H",
    ("brca1_bard1", "K45R"): "Y2H",
    ("mlh1_pms2", "L749P"): "co-IP",
    ("mlh1_pms2", "Q542L"): "co-IP",
    ("mlh1_pms2", "R755W"): "co-IP; MAVE interaction score",
    ("msh2_msh6", "N127S"): "MSH6 IP",
    ("smad4_smad3", "D351H"): "co-IP",
    ("smad4_smad3", "R361C"): "co-IP",
    ("smad4_smad3", "I500V"): "GST pull-down; TR-FRET",
}
# Tripwire: the variants the rule reaches. A ledger or canonical change that
# alters this set must be reviewed, not absorbed.
EXPECTED_EITHER_AXIS = {
    ("brca1_bard1", "V11G"), ("mlh1_pms2", "L749P"), ("mlh1_pms2", "Q542L"),
    ("smad4_smad3", "D351H"), ("smad4_smad3", "R361C"), ("smad4_smad3", "I500V"),
}
ASSAY_PATTERN = re.compile(
    r"co-?IP|immunoprecip|\bIP\b|pull-?down|two-hybrid|\bY2H\b|TR-FRET", re.I)
DIRECTIONAL = {"destab", "stab"}


def load_ledger() -> pd.DataFrame:
    return pd.read_csv(LEDGER, keep_default_na=False, dtype=str)


def curated_ledger(ledger: pd.DataFrame) -> pd.DataFrame:
    """The ledger with complex_amount_assay set from COMPLEX_AMOUNT."""
    out = ledger.copy()
    out[LEDGER_COLUMN] = ""
    seen = set()
    for idx, row in out.iterrows():
        key = (row.system, row.variant)
        text = "%s %s" % (row.evidence_basis, row.evidence_citation)
        if row.axis == "binding" and key in COMPLEX_AMOUNT:
            for name in COMPLEX_AMOUNT[key].split("; "):
                if name.lower() not in text.lower():
                    raise SystemExit("STOP: %s %s: %r not named in the basis or citation"
                                     % (key[0], key[1], name))
            out.at[idx, LEDGER_COLUMN] = COMPLEX_AMOUNT[key]
            seen.add(key)
        elif row.axis == "binding" and ASSAY_PATTERN.search(text):
            raise SystemExit("STOP: %s %s binding names an interaction assay but is not "
                             "classified in COMPLEX_AMOUNT: %s" % (key[0], key[1], text[:120]))
    missing = set(COMPLEX_AMOUNT) - seen
    if missing:
        raise SystemExit("STOP: COMPLEX_AMOUNT entries with no binding row: %s" % sorted(missing))
    return out


def either_axis_variants(ledger: pd.DataFrame, canon: pd.DataFrame) -> set:
    """Variants whose directional binding expectation rests on a complex-amount
    assay and whose assembled-complex stability axis was not tested."""
    cmap = canon.set_index(["system", "variant"])
    found = set()
    for row in ledger.itertuples():
        if row.axis != "binding" or not getattr(row, LEDGER_COLUMN):
            continue
        if row.expected_token not in DIRECTIONAL:
            continue
        key = (row.system, row.variant)
        crow = cmap.loc[key]
        if str(crow.expected_ddg_binding).strip().lower() != row.expected_token:
            raise SystemExit("STOP: %s %s binding token disagrees between ledger and canonical" % key)
        # The rule needs an assembled-complex stability axis that exists and was
        # not tested ("unknown"); a committed or not_applicable token is reviewed.
        fold_c = str(crow.expected_ddg_fold_complex).strip().lower()
        if fold_c != "unknown":
            raise SystemExit("STOP: %s %s assembled-complex token is %r, not untested; "
                             "review before applying the rule" % (key[0], key[1], fold_c))
        if crow.expected_mech_class == "structurally_silent":
            raise SystemExit("STOP: %s %s is non-disruptive but carries a directional "
                             "binding expectation" % key)
        found.add(key)
    if found != EXPECTED_EITHER_AXIS:
        raise SystemExit("STOP: either_axis set changed: +%s -%s"
                         % (sorted(found - EXPECTED_EITHER_AXIS), sorted(EXPECTED_EITHER_AXIS - found)))
    return found


def with_scope(canon: pd.DataFrame, either: set) -> pd.DataFrame:
    """The canonical with the scope column set, placed immediately after
    expected_topology so it sits with the ground truth it qualifies. A fixed
    position also keeps the canonical updater idempotent: it re-appends its
    ISDS columns at the end, so a column appended after them would move.
    """
    out = canon.drop(columns=[SCOPE_COLUMN], errors="ignore").copy()
    values = [ac.EITHER_ASSEMBLED_AXIS if (s, v) in either else ""
              for s, v in zip(out.system, out.variant)]
    out.insert(out.columns.get_loc("expected_topology") + 1, SCOPE_COLUMN, values)
    return out


def _same(a, b) -> bool:
    if pd.isna(a) and pd.isna(b):
        return True
    if pd.isna(a) or pd.isna(b):
        return False
    return abs(float(a) - float(b)) < 1e-12


def agreement(canon: pd.DataFrame, either: set):
    """Recompute the stored per-variant agreement columns at every threshold.

    Per-axis direction agreement scores a complex-amount expectation jointly on
    either assembled axis (apply_concordance_v5.structural_agreement_by_axis,
    'assembled_either' row; compute_directional_agreement likewise). The stored
    columns were last written in v7.6, so they are recomputed here with the
    pipeline's own functions, in the pipeline's own formulas. Any change outside
    the either_axis variants stops the run: it would mean the stored columns
    were stale for another reason, which must be investigated, not absorbed.
    """
    out = canon.copy()
    partners = ac.discover_partners(canon)
    changes = []
    for tag, spec in ac.THRESHOLD_SPECS:
        mt, ft, bt = ((spec["monomer"], spec["fold"], spec["binding"])
                      if isinstance(spec, dict) else (float(spec),) * 3)
        for idx, row in canon.iterrows():
            sn, sd = ac.compute_structural_agreement(row, partners, mt, ft, bt)
            df_, dh, dd = ac.compute_directional_agreement(row, partners, mt, ft, bt)
            new = {
                "structural_agreement_n_%s" % tag: sn,
                "structural_agreement_d_%s" % tag: sd,
                "structural_agreement_%s" % tag: (sn / sd if sd > 0 else None),
                "directional_agreement_full_%s" % tag: df_,
                "directional_agreement_half_%s" % tag: dh,
                "directional_agreement_d_%s" % tag: dd,
                "directional_agreement_%s" % tag: ((df_ + 0.5 * dh) / dd if dd > 0 else None),
            }
            for col, value in new.items():
                if col not in canon.columns:
                    raise SystemExit("STOP: canonical lacks %s" % col)
                old = row[col]
                if not _same(old, value):
                    changes.append({"system": row.system, "variant": row.variant,
                                    "column": col, "old": old, "new": value})
                    out.at[idx, col] = value
    changes = pd.DataFrame(changes, columns=["system", "variant", "column", "old", "new"])
    stray = changes[[(s, v) not in either for s, v in zip(changes.system, changes.variant)]]
    if len(stray):
        raise SystemExit("STOP: agreement columns changed outside the either_axis variants:\n%s"
                         % stray.head(12).to_string())
    return out, changes


def annotate(changes: pd.DataFrame, either: set) -> pd.DataFrame:
    if not len(changes):
        return changes
    stray = changes[[(s, v) not in either for s, v in zip(changes.system, changes.variant)]]
    if len(stray):
        raise SystemExit("STOP: grades changed outside the either_axis variants:\n%s"
                         % stray[["pipeline", "threshold", "variant", "old_grade", "new_grade"]].to_string())
    out = changes.copy()
    out["assay"] = [COMPLEX_AMOUNT[(s, v)] for s, v in zip(out.system, out.variant)]
    out["reason"] = [
        "evidence resolution: the %s measures how much complex forms, so a call in the "
        "measured direction on either assembled axis matches it" % a for a in out.assay]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    ledger = load_ledger()
    curated = curated_ledger(ledger.drop(columns=[LEDGER_COLUMN], errors="ignore"))
    canon = pd.read_csv(CANON, low_memory=False)
    either = either_axis_variants(curated, canon)
    scoped = with_scope(canon, either)
    regraded, changes = table3.regrade(scoped)
    changes = annotate(changes, either)
    agreed, agree_changes = agreement(regraded, either)

    if args.check:
        failures = []
        if LEDGER_COLUMN not in ledger.columns or not ledger[LEDGER_COLUMN].equals(curated[LEDGER_COLUMN]):
            failures.append("ledger %s does not match the curation" % LEDGER_COLUMN)
        stored = canon[SCOPE_COLUMN].fillna("").astype(str) if SCOPE_COLUMN in canon.columns else None
        if stored is None or list(stored) != list(scoped[SCOPE_COLUMN]):
            failures.append("canonical %s does not match the ledger rule" % SCOPE_COLUMN)
        elif list(canon.columns) != list(scoped.columns):
            failures.append("canonical %s is not placed after expected_topology" % SCOPE_COLUMN)
        if len(changes):
            failures.append("%d stored grades disagree with the evidence-resolution grader" % len(changes))
        if len(agree_changes):
            failures.append("%d stored agreement cells disagree with the joint scoring"
                            % len(agree_changes))
        for f in failures:
            print("FAIL: %s" % f)
        if failures:
            return 1
        print("PASS: ledger column, canonical scope (%d variants), every stored grade and every "
              "stored agreement cell match the evidence-resolution rule" % len(either))
        return 0

    applied = (LEDGER_COLUMN in ledger.columns and ledger[LEDGER_COLUMN].equals(curated[LEDGER_COLUMN])
               and SCOPE_COLUMN in canon.columns
               and list(canon[SCOPE_COLUMN].fillna("").astype(str)) == list(scoped[SCOPE_COLUMN])
               and list(canon.columns) == list(scoped.columns))
    if applied and not len(changes) and not len(agree_changes):
        print("already applied; nothing to do (the audit files are left as written)")
        return 0
    touched = [c for c in canon.columns if not canon[c].equals(agreed[c])]
    allowed = ((SCOPE_COLUMN, "structural_agreement_", "directional_agreement_")
               + tuple("%smech_%s" % (p, s) for p, _ in table3.PIPELINES
                       for s in ("consistency_", "false_positive_axes_", "missed_positive_axes_")))
    stray = [c for c in touched if not c.startswith(allowed)]
    if stray:
        raise SystemExit("STOP: columns outside the grade and agreement families would change: %s" % stray)
    curated.to_csv(LEDGER, index=False)
    agreed.to_csv(CANON, index=False)
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    # Each audit is written only when its part of the rule changed something,
    # so re-running after a partial application never erases an earlier audit.
    if len(changes):
        changes.to_csv(AUDIT, index=False)
    if len(agree_changes):
        agree_changes.to_csv(AGREEMENT_AUDIT, index=False)
    print("ledger: %s set on %d binding rows" % (LEDGER_COLUMN, int((curated[LEDGER_COLUMN] != "").sum())))
    print("canonical: %s = %s for %d variants: %s" % (SCOPE_COLUMN, ac.EITHER_ASSEMBLED_AXIS, len(either),
                                                     ", ".join("%s %s" % k for k in sorted(either))))
    print("grades changed: %d%s" % (len(changes), " (written to %s)" % AUDIT.relative_to(REPO)
                                     if len(changes) else ""))
    if len(changes):
        print(changes.groupby(["pipeline", "threshold"]).size().to_string())
    print("agreement cells changed: %d%s" % (len(agree_changes), " (written to %s)"
                                             % AGREEMENT_AUDIT.relative_to(REPO) if len(agree_changes) else ""))
    if len(agree_changes):
        print(agree_changes.groupby("variant").size().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
