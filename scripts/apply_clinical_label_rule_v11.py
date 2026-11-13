#!/usr/bin/env python3
"""Clinical labels follow the Franklin classification, applied uniformly.

WHY THIS EXISTS
---------------
Table 1 says its Pathogenic column counts the variants "annotated pathogenic in
the clinical databases", and the clinical labels feed every pathogenicity
comparison in the paper. Before v11, 38 of the 46 labelled benchmark variants
in complexes followed their Franklin classification, and 8 followed neither
Franklin nor ClinVar:

  VWF R1334Q        pathogenic   Franklin VUS        ClinVar uncertain significance
  HBB W37Y/A/E      pathogenic   Franklin none       no ClinVar record
  HBB W37G          pathogenic   Franklin none       ClinVar uncertain significance
  HBB N102T         pathogenic   Franklin VUS-low    ClinVar benign (no assertion criteria)
  BRCA1 E100Q       benign       Franklin VUS-mid    ClinVar uncertain significance
  PMS2 G857A        benign       Franklin VUS-mid    ClinVar benign (expert panel)

From v11 the clinical label follows one rule, in order of precedence:
  1. A ClinVar classification reviewed by an expert panel or given in a
     practice guideline (data/analysis_inputs/clinvar/COMAVI_clinvar_records.csv,
     retrieved by scripts/fetch_clinvar_records.py). If expert-panel records for
     the same protein change disagree across nucleotide changes, this step is
     skipped (PIK3R1 M326I: c.978G>A Benign, c.978G>C Uncertain significance).
  2. Otherwise the Franklin classification.
Pathogenic or likely pathogenic -> pathogenic (gain of function kept where the
benchmark records it); benign or likely benign -> benign; uncertain significance
or any VUS grade -> vus; no classification -> unclassified. Only pathogenic and
benign variants enter a clinical comparison.

Expert-panel precedence changes one label relative to Franklin alone: PMS2
G857A, Franklin VUS-mid, is Benign by a ClinVar expert panel (VCV000036691),
which restores its v10.1 label. Twelve other variants carry an expert-panel
classification and every one agrees with Franklin; PIK3R1 M326I carries two that
disagree with each other (step 1 above), so Franklin decides it.

A Franklin classification that contradicts an expert panel, or a rule label
that reverses a v10.1 label (pathogenic vs benign), stops the run, because that
would be a curation question, not a rule.

The audit (docs/audit/COMAVI_v11_clinical_label_rule.csv) is rebuilt in full on
every run: every variant whose label differs from v10.1, and every variant whose
label rests on an expert-panel record.

Usage:
    python scripts/apply_clinical_label_rule_v11.py          # apply, write the audit
    python scripts/apply_clinical_label_rule_v11.py --check  # verify only
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
AUDIT = REPO / "docs" / "audit" / "COMAVI_v11_clinical_label_rule.csv"

CLINVAR = REPO / "data" / "analysis_inputs" / "clinvar" / "COMAVI_clinvar_records.csv"
EXPERT_REVIEW = ("reviewed by expert panel", "practice guideline")

LABELLED = ("pathogenic", "pathogenic_gof", "benign")
# The v10.1 labels of the variants whose label the rule first changed (commit
# e1dcbb3, reference_outputs/scored_61var_canonical.csv, column phenotype).
# Every other variant kept its v10.1 label. Frozen here so the audit can be
# rebuilt in full on every run rather than only when a label moves.
PRE_V11_LABELS = {
    ("vwf_gpiba", "R1334Q"): "pathogenic",
    ("hemoglobin_tetramer", "W37Y"): "pathogenic",
    ("hemoglobin_tetramer", "W37A"): "pathogenic",
    ("hemoglobin_tetramer", "W37E"): "pathogenic",
    ("hemoglobin_tetramer", "W37G"): "pathogenic",
    ("hemoglobin_tetramer", "N102T"): "pathogenic",
    ("brca1_bard1", "E100Q"): "benign",
    ("mlh1_pms2", "G857A"): "benign",
}


def franklin_class(value) -> str:
    v = str(value).strip().lower() if pd.notna(value) else ""
    if v in ("pathogenic", "likely_pathogenic"):
        return "pathogenic"
    if v in ("benign", "likely_benign"):
        return "benign"
    if v.startswith("vus"):
        return "vus"
    return "unclassified"


def clinvar_class(description) -> str:
    v = str(description).strip().lower()
    if v in ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"):
        return "pathogenic"
    if v in ("benign", "likely benign", "benign/likely benign"):
        return "benign"
    if v == "uncertain significance":
        return "vus"
    raise SystemExit("STOP: unmapped ClinVar classification %r" % description)


def expert_panel_classes() -> dict:
    """{(system, variant): (class or None, note)} from the committed ClinVar table.

    Only records reviewed by an expert panel or given in a practice guideline
    count. When such records for the same protein change disagree (different
    nucleotide changes can encode one substitution), the rule falls back to
    Franklin and the note says why.
    """
    table = pd.read_csv(CLINVAR, keep_default_na=False, dtype=str)
    expert = table[table.review_status.isin(EXPERT_REVIEW)]
    out = {}
    for (system, variant), grp in expert.groupby(["system", "variant"]):
        classes = {clinvar_class(c) for c in grp.classification}
        vcvs = ", ".join(sorted(grp.vcv))
        if len(classes) == 1:
            out[(system, variant)] = (classes.pop(), "ClinVar expert panel (%s)" % vcvs)
        else:
            out[(system, variant)] = (None, "expert-panel records disagree across nucleotide "
                                            "changes (%s); Franklin used" % vcvs)
    return out


def rule_label(stored, franklin, expert=None) -> str:
    f = franklin_class(franklin)
    if expert is not None:
        if {f, expert} == {"pathogenic", "benign"}:
            raise SystemExit("STOP: Franklin %s contradicts a ClinVar expert panel (%s) -- a "
                             "curation question, not a rule" % (franklin, expert))
        f = expert
    if f == "pathogenic" and stored == "pathogenic_gof":
        return "pathogenic_gof"
    return f


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    canon = pd.read_csv(CANON, low_memory=False)
    panels = expert_panel_classes()
    rows, bad = [], []
    for i, r in canon.iterrows():
        stored = r.phenotype
        if stored not in LABELLED and stored not in ("vus", "unclassified"):
            continue                          # BRCT and other unannotated rows
        key = (r.system, r.variant)
        expert, note = panels.get(key, (None, ""))
        new = rule_label(stored, r.franklin, expert)
        before = PRE_V11_LABELS.get(key, stored)
        if before in LABELLED and new in LABELLED and (before == "benign") != (new == "benign"):
            raise SystemExit("STOP: %s %s was %s before v11 but the rule gives %s -- a curation "
                             "question, not a rule" % (r.system, r.variant, before, new))
        source = note if expert is not None else ("Franklin" + ("; " + note if note else ""))
        if new != before or expert is not None and franklin_class(r.franklin) != expert or note:
            rows.append(dict(system=r.system, gene=r.gene, variant=r.variant,
                             label_before_v11=before, franklin=r.franklin, label_source=source,
                             label_v11=new))
        if new != stored:
            if args.check:
                bad.append("%s %s: stored %s, rule gives %s" % (r.system, r.variant, stored, new))
            else:
                canon.at[i, "phenotype"] = new
    audit = pd.DataFrame(rows)
    if args.check:
        if not AUDIT.exists() or not pd.read_csv(AUDIT, keep_default_na=False).astype(str).equals(
                audit.fillna("").astype(str)):
            bad.append("the audit %s is not current" % AUDIT.name)
        for b in bad:
            print("  MISMATCH " + b)
        print("%s: clinical labels %s the rule (expert panel, then Franklin)"
              % ("FAIL" if bad else "PASS", "disagree with" if bad else "follow"))
        return 1 if bad else 0
    canon.to_csv(CANON, index=False)
    audit.to_csv(AUDIT, index=False)
    print("audit rows (every label that differs from v10.1 or rests on ClinVar): %d" % len(audit))
    print(audit[["gene", "variant", "label_before_v11", "franklin", "label_v11", "label_source"]]
          .to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
