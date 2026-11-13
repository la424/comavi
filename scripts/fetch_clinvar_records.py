#!/usr/bin/env python3
"""Retrieve the ClinVar records for the benchmark variants in complexes.

WHY THIS EXISTS
---------------
From v11 the clinical label of a benchmark variant follows a fixed rule
(scripts/apply_clinical_label_rule_v11.py): a ClinVar classification reviewed by
an expert panel or given in a practice guideline takes precedence, and otherwise
the Franklin classification is used. This script retrieves the ClinVar records
the rule reads, once, through NCBI E-utilities, and writes them to a committed
table so the labels can be re-derived offline and audited record by record.
ClinVar changes over time; the retrieval date is recorded with every row, and
re-running this script is a deliberate update, not part of CI.

QUERY
-----
For each variant: esearch db=clinvar, term '<GENE>[gene] AND "p.<Ref><pos><Alt>"'
(three-letter HGVS), then esummary, keeping records whose title names exactly
that protein change. Numbering rules, as in scripts/apply_alphamissense_release.py:
  * HBB is numbered from the mature chain in the benchmark, so HGVS = +1.
  * CFH I62V and VWF A1381T are named in the reverse direction in the benchmark
    (the reference residue is the second letter).
  * Calmodulin is encoded by three genes with an identical protein. When the
    benchmark gene has no record for the change, CALM2 and CALM3 are searched
    and the gene that holds the record is written to the row (CALM1 D96V, the
    variant of Crotti 2013, is recorded under CALM2).

Usage:
    python scripts/fetch_clinvar_records.py      # writes the table (network)
"""
from __future__ import annotations

import datetime
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import benchmark_population as bp  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
OUT = REPO / "data" / "analysis_inputs" / "clinvar" / "COMAVI_clinvar_records.csv"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
       "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
       "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}
MATURE_NUMBERING = {"HBB"}
REVERSE_NAMED = {("CFH", "I62V"), ("VWF", "A1381T")}
PARALOGS = {"CALM1": ("CALM2", "CALM3")}


def hgvs_protein(gene: str, variant: str) -> str:
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", variant)
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    if (gene, variant) in REVERSE_NAMED:
        ref, alt = alt, ref
    if gene in MATURE_NUMBERING:
        pos += 1
    return "p.%s%d%s" % (AA3[ref], pos, AA3[alt])


def get_json(url: str):
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=40) as resp:
                return json.load(resp)
        except Exception as exc:  # transient network or rate-limit errors
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


def records(gene: str, hgvs_p: str) -> list[dict]:
    term = urllib.parse.quote('%s[gene] AND "%s"' % (gene, hgvs_p))
    ids = get_json("%s/esearch.fcgi?db=clinvar&term=%s&retmode=json&retmax=40"
                   % (EUTILS, term))["esearchresult"]["idlist"]
    time.sleep(0.4)
    if not ids:
        return []
    res = get_json("%s/esummary.fcgi?db=clinvar&id=%s&retmode=json" % (EUTILS, ",".join(ids)))["result"]
    time.sleep(0.4)
    out = []
    for uid in res["uids"]:
        rec = res[uid]
        title = rec.get("title", "")
        if "(%s)" % hgvs_p not in title:
            continue
        cls = rec.get("germline_classification") or {}
        out.append({
            "clinvar_gene": gene, "vcv": rec.get("accession", ""), "title": title,
            "classification": cls.get("description", ""), "review_status": cls.get("review_status", ""),
            "last_evaluated": cls.get("last_evaluated", ""),
            "conditions": "; ".join(t.get("trait_name", "") for t in (cls.get("trait_set") or []))})
    return out


def main() -> int:
    canon = pd.read_csv(CANON, low_memory=False)
    bench = bp.benchmark(canon)
    complexes = bench[bench.system != "brca1_brct"]
    retrieved = datetime.date.today().isoformat()
    rows = []
    for r in complexes.itertuples():
        gene = str(r.gene).upper()
        hgvs_p = hgvs_protein(gene, r.variant)
        found = records(gene, hgvs_p)
        if not found:
            for alt in PARALOGS.get(gene, ()):
                found = records(alt, hgvs_p)
                if found:
                    break
        base = {"system": r.system, "gene": gene, "variant": r.variant, "hgvs_p": hgvs_p,
                "retrieved": retrieved}
        if not found:
            rows.append({**base, "clinvar_gene": "", "vcv": "", "title": "", "classification": "",
                         "review_status": "no record", "last_evaluated": "", "conditions": ""})
        for rec in found:
            rows.append({**base, **rec})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    print("wrote %s: %d variants, %d records, %d without a record"
          % (OUT.relative_to(REPO), table[["system", "variant"]].drop_duplicates().shape[0],
             int((table.review_status != "no record").sum()), int((table.review_status == "no record").sum())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
