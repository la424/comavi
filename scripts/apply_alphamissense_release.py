#!/usr/bin/env python3
"""Set every AlphaMissense value in the benchmark from the checksummed release.

WHY THIS EXISTS
---------------
Before v11.0 the AlphaMissense column had two undocumented origins: 44 values
typed into inputs/AM_variants_comavi_mechanism_test.xlsx at three decimals, and
three CFH values written straight into the canonical table when the cohort was
expanded. The two VWF variants had none, which is why the pathogenicity
comparison could not use the prioritization population.

Checked against the release (scripts/fetch_alphamissense_subset.py), 46 of the
47 existing values reproduce at three decimals. The forty-seventh, MLH1 V384D,
read 0.988 where the release gives 0.9976; 0.988 is the release value for the
neighbouring substitution V384E, so the value was transcribed from the wrong
row. V384D is clinically benign, and the error raised the published
AlphaMissense pathogenicity AUC from 0.886 to 0.903 on the 47-variant set.

From v11.0 every value comes from one file whose MD5 matches the one published
with the record, at its full four-decimal precision, by an explicit rule:

  * UniProt canonical numbering from the initiator methionine, except HBB,
    whose variants follow the globin convention that excludes it (+1);
  * two calculations are named in the reverse direction because the deposited
    structure carries the substituted residue (CFH I62V is clinical
    p.Val62Ile; VWF A1381T is clinical p.Thr1381Ala), so the release row is
    the clinical substitution.

Every lookup also confirms the release row's reference residue, which the
release key encodes.

Usage:
    python scripts/apply_alphamissense_release.py          # apply and write the audit
    python scripts/apply_alphamissense_release.py --check  # verify, change nothing
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import openpyxl
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
SUBSET = REPO / "inputs" / "alphamissense" / "AlphaMissense_benchmark_subset.tsv"
PROV = REPO / "inputs" / "alphamissense" / "AlphaMissense_benchmark_subset.provenance.json"
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
XLSX = REPO / "inputs" / "AM_variants_comavi_mechanism_test.xlsx"
AUDIT = REPO / "docs" / "audit" / "COMAVI_alphamissense_source_reconciliation.csv"

# Every variant in a complex that carries or carried a clinical annotation.
# v11: "vus" and "unclassified" (scripts/apply_clinical_label_rule_v11.py)
# keep their scores; they no longer enter a clinical comparison.
CLINICAL = ("pathogenic", "pathogenic_gof", "benign", "vus", "unclassified")
CLASS_NAMES = {"pathogenic": "likely_pathogenic", "benign": "likely_benign",
               "ambiguous": "ambiguous"}
MATURE_NUMBERING = {"HBB"}
REVERSE_NAMED = {("CFH", "I62V"), ("VWF", "A1381T")}
AA3 = dict(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu", G="Gly",
           H="His", I="Ile", L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser",
           T="Thr", W="Trp", Y="Tyr", V="Val")


def load_release():
    prov = json.loads(PROV.read_text())
    if not prov.get("md5_verified"):
        raise SystemExit("release subset was not checksum-verified; re-run "
                         "scripts/fetch_alphamissense_subset.py")
    sub = pd.read_csv(SUBSET, sep="\t")
    lut = {(a, v): (float(s), c) for a, v, s, c in
           zip(sub.uniprot_id, sub.protein_variant, sub.am_pathogenicity, sub.am_class)}
    return prov, lut


def release_key(gene: str, variant: str):
    g = gene.upper()
    ref, pos, alt = variant[0], int(variant[1:-1]), variant[-1]
    rule = "as named"
    if g in MATURE_NUMBERING:
        pos += 1
        rule = "globin numbering (+1 for the initiator Met)"
    if (g, variant) in REVERSE_NAMED:
        ref, alt = alt, ref
        rule = "reverse-named calculation"
    return "%s%d%s" % (ref, pos, alt), rule


def lookup(prov, lut, gene, variant):
    acc = prov["accessions"][gene.upper()]
    key, rule = release_key(gene, variant)
    if (acc, key) not in lut:
        raise KeyError("%s %s -> %s %s is not in the release (wrong reference "
                       "residue or numbering)" % (gene, variant, acc, key))
    score, cls = lut[(acc, key)]
    hgvs = "p.%s%s%s" % (AA3[key[0]], key[1:-1], AA3[key[-1]])
    return acc, key, rule, hgvs, round(score, 4), CLASS_NAMES[cls]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    prov, lut = load_release()

    canon = pd.read_csv(CANON, low_memory=False)
    rows, bad = [], []
    for i, r in canon.iterrows():
        if r.phenotype not in CLINICAL:
            continue
        acc, key, rule, hgvs, score, cls = lookup(prov, lut, str(r.gene), r.variant)
        prev, prev_cls = r["AM pathogenicity"], r["AM class"]
        if pd.isna(prev):
            status = "added"
        elif abs(round(float(prev), 3) - round(score, 3)) <= 0.0011:
            status = "unchanged at three decimals"
        else:
            status = "corrected"
        rows.append(dict(system=r.system, gene=str(r.gene).upper(), variant=r.variant,
                         clinical_substitution=hgvs, uniprot=acc, release_key=key,
                         numbering_rule=rule, previous_value=prev, release_value=score,
                         previous_class=prev_cls, release_class=cls, status=status))
        if args.check:
            if pd.isna(prev) or abs(float(prev) - score) > 1e-9 or prev_cls != cls:
                bad.append("%s %s: canonical %s/%s, release %s/%s"
                           % (r.system, r.variant, prev, prev_cls, score, cls))
        else:
            canon.at[i, "AM pathogenicity"] = score
            canon.at[i, "AM class"] = cls

    wb = openpyxl.load_workbook(XLSX)
    ws = wb.active
    head = [c.value for c in ws[1]]
    gi, vi = head.index("gene"), head.index("variant")
    si, ci = head.index("AM pathogenicity"), head.index("AM class")
    for row in ws.iter_rows(min_row=2):
        gene, variant = row[gi].value, row[vi].value
        if not gene or not variant:
            continue
        _, _, _, _, score, cls = lookup(prov, lut, str(gene), str(variant))
        if args.check:
            if abs(float(row[si].value) - score) > 1e-9 or row[ci].value != cls:
                bad.append("xlsx %s %s: %s/%s, release %s/%s"
                           % (gene, variant, row[si].value, row[ci].value, score, cls))
        else:
            row[si].value, row[ci].value = score, cls

    if args.check:
        for b in bad:
            print("  MISMATCH " + b)
        print("%s: %d canonical rows and %d input rows checked against the release"
              % ("FAIL" if bad else "PASS", len(rows), ws.max_row - 1))
        return 1 if bad else 0

    canon.to_csv(CANON, index=False)
    wb.save(XLSX)
    audit = pd.DataFrame(rows)
    audit.to_csv(AUDIT, index=False)
    print("applied: %s" % audit.status.value_counts().to_dict())
    for _, a in audit[audit.status != "unchanged at three decimals"].iterrows():
        print("  %-10s %-7s %-14s %s -> %s (%s)" % (a.system, a.variant, a.clinical_substitution,
                                                  a.previous_value, a.release_value, a.status))
    return 0


if __name__ == "__main__":
    sys.exit(main())
