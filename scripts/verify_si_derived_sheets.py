#!/usr/bin/env python3
"""Recompute the supplementary sheets that hold derived values.

WHY THIS EXISTS
---------------
Nothing checked the contents of the supplementary workbooks. The final audit
demonstrated this twice over.

First, by finding a real defect the hard way: the S3 Table tier_distribution
sheet still read "Tier 2  6/17 = 0.353" and asserted "Tiers 3-4  0/17 = 0.000"
long after the ground truth moved, while the S2 Text prose next to it already
said 4/17. The supplement contradicted itself on a perfect-value claim, and the
only thing that caught it was enumerating every cell by hand.

Second, by negative control: editing a cell of S3_Table.xlsx to read
"11/13 = 0.999" and re-running verify_si_package.py, audit_evidence_claims.py
and verify_denominators.py produced three passes. A workbook cell was, until
this file, unreachable by any gate in the repository.

A checksum manifest is not a substitute. The stale sheet was already stale when
the package was built, so its digest matched perfectly. The only check that
finds this class is recomputing the cell from the record it is derived from.

Each sheet below names the record it is recomputed from. Sheets that are
transcriptions of a record (threshold_sweep, robustness) are compared value by
value; sheets that are curated prose are out of scope and listed at the end.

Usage:
    PYTHONPATH=. python scripts/verify_si_derived_sheets.py
"""
from __future__ import annotations

import glob
import json
import pathlib
import re
import sys

import openpyxl
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"


def si_dir() -> pathlib.Path:
    """Resolve the SI directory by glob, never by a versioned literal.

    A hardcoded directory name is how an earlier extraction script came to
    contribute nothing while still reporting success: the name moved and its
    loop simply matched no files.
    """
    hits = sorted(glob.glob(str(REPO / "submission" / "COMAVI_*_SI")))
    if len(hits) != 1:
        raise SystemExit("FAIL: expected exactly one submission/COMAVI_*_SI "
                         "directory, found %d: %s" % (len(hits), hits))
    return pathlib.Path(hits[0])


def rows_of(path: pathlib.Path, sheet: str) -> list[list[str]]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        if sheet not in wb.sheetnames:
            raise SystemExit("FAIL: sheet %s absent from %s" % (sheet, path.name))
        out = []
        for r in wb[sheet].iter_rows(values_only=True):
            if r and any(c is not None for c in r):
                out.append(["" if c is None else str(c).strip() for c in r])
        return out
    finally:
        wb.close()


def num(s: str):
    m = re.search(r"-?\d+(?:\.\d+)?", str(s))
    return float(m.group(0)) if m else None


FAILURES: list[str] = []
CHECKED: list[str] = []


def main() -> int:
    SI = si_dir()
    s3 = SI / "S3_Table.xlsx"

    # ---- tier_distribution, from ISDS_v1_per_variant.csv
    pv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
    tn = pv.reported_tier.astype(str).str.extract(r"(\d)").astype(float)[0]
    y = pv.structural_ground_truth.astype(bool)
    want = {}
    for label, mask in (("Tier 1", tn == 1), ("Tier 2", tn == 2),
                        ("Tiers 3-4", tn.isin([3, 4]))):
        k = int((mask & y).sum()); n = int(mask.sum())
        want[label] = (k, n - k, "%d/%d = %.3f" % (k, n, k / n if n else 0.0))
    for row in rows_of(s3, "tier_distribution"):
        if row and row[0] in want:
            k, nonk, frac = want[row[0]]
            got = (num(row[1]), num(row[2]), row[3])
            if (int(got[0]), int(got[1])) != (k, nonk) or got[2] != frac:
                FAILURES.append(
                    "tier_distribution %s: sheet says %s / %s / %s, record gives "
                    "%d / %d / %s" % (row[0], row[1], row[2], row[3], k, nonk, frac))
    CHECKED.append("tier_distribution <- ISDS_v1_per_variant.csv")

    # ---- threshold_sweep, from COMAVI_threshold_sweep.json
    sweep = {r["threshold_kcal_mol"] if "threshold_kcal_mol" in r else r["threshold_tag"]: r
             for r in json.loads((RO / "COMAVI_threshold_sweep.json").read_text())["rows"]}
    by_kcal = {}
    for r in json.loads((RO / "COMAVI_threshold_sweep.json").read_text())["rows"]:
        key = str(r.get("threshold_kcal_mol", r["threshold_tag"]))
        by_kcal[key] = r
    for row in rows_of(s3, "threshold_sweep"):
        rec = by_kcal.get(row[0])
        if not rec:
            continue
        pairs = [("Detection", "detection"), ("Attribution", "attribution"),
                 ("Correct rejection", "correct_rejection"),
                 ("Whole-variant score", "whole_variant")]
        for col, (_, field) in zip(row[1:5], pairs):
            got, exp = num(col), rec.get(field)
            if got is None or exp is None:
                continue
            if abs(got - float(exp)) > 5e-4:
                FAILURES.append("threshold_sweep %s %s: sheet %.4f, record %.4f"
                                % (row[0], field, got, float(exp)))
    CHECKED.append("threshold_sweep <- COMAVI_threshold_sweep.json")

    # ---- robustness: the observed row, from the sweep and the numbers ledger
    nl = json.loads((RO / "COMAVI_numbers_ledger.json").read_text())
    t25 = by_kcal.get("2.5") or by_kcal.get("t25")
    for row in rows_of(s3, "robustness"):
        if row and row[0].lower().startswith("observed"):
            got_mc, got_sa = num(row[1]), num(row[2])
            exp_mc = float(t25["whole_variant"])
            exp_sa = nl.get("SA_total")
            if got_mc is not None and abs(got_mc - exp_mc) > 5e-4:
                FAILURES.append("robustness observed MC: sheet %.4f, record %.4f"
                                % (got_mc, exp_mc))
            if got_sa is not None and isinstance(exp_sa, (int, float)) \
                    and abs(got_sa - float(exp_sa)) > 5e-4:
                FAILURES.append("robustness observed SA: sheet %.4f, record %.4f"
                                % (got_sa, float(exp_sa)))
    CHECKED.append("robustness (observed row) <- threshold sweep + numbers ledger")

    # ---- S1 evidence_ledger, from COMAVI_evidence_ledger.csv
    #
    # This sheet was classified OUT OF SCOPE on the first pass, as curated
    # rather than derived. That was wrong: it is a COPY of the committed
    # ledger, so it is exactly as derivable as any other sheet here, and
    # leaving it ungated let it drift through two correction rounds. When the
    # gap was finally found the shipped sheet was still in its pre-v8.0 state:
    # it was missing both rows v8.0 added (mlh1_pms2 H718Y monomer,
    # smad4_smad3 I500V binding), still carried the row v8.0 retired
    # (msh2_msh6 G674R binding), and five evidence_basis cells held superseded
    # text. The supplement therefore advertised 98 committed axes against the
    # article's stated 99, and published a commitment the correction had
    # explicitly withdrawn. Nothing caught it, because every other gate reads
    # records rather than the workbook.
    led = pd.read_csv(RO / "COMAVI_evidence_ledger.csv")
    s1_rows = rows_of(SI / "S1_Table.xlsx", "evidence_ledger")
    head = next((i for i, r in enumerate(s1_rows)
                 if r[:3] == ["system", "variant", "axis"]), None)
    if head is None:
        FAILURES.append("S1 evidence_ledger: no header row found")
    else:
        cols = s1_rows[head]
        body = [r for r in s1_rows[head + 1:] if r and r[0]]
        sheet_keys = {(r[0], r[1], r[2]) for r in body}
        led_keys = set(zip(led.system, led.variant, led.axis))
        for k in sorted(led_keys - sheet_keys):
            FAILURES.append("S1 evidence_ledger: committed axis missing from "
                            "the sheet: %s %s %s" % k)
        for k in sorted(sheet_keys - led_keys):
            FAILURES.append("S1 evidence_ledger: sheet carries an axis the "
                            "ledger does not commit: %s %s %s" % k)
        bi = cols.index("evidence_basis")
        ti = cols.index("expected_token")
        auth = {(r.system, r.variant, r.axis): (str(r.evidence_basis),
                                                str(r.expected_token))
                for _, r in led.iterrows()}
        for r in body:
            k = (r[0], r[1], r[2])
            if k not in auth:
                continue
            basis, token = auth[k]
            if r[bi].strip() != basis.strip():
                FAILURES.append("S1 evidence_ledger %s %s %s: basis text "
                                "differs from the ledger" % k)
            if r[ti].strip() != token.strip():
                FAILURES.append("S1 evidence_ledger %s %s %s: token '%s' "
                                "against ledger '%s'" % (k + (r[ti], token)))
    CHECKED.append("S1 evidence_ledger <- COMAVI_evidence_ledger.csv "
                   "(membership, token and basis text)")

    if FAILURES:
        print("FAIL: %d supplementary cells disagree with their record" % len(FAILURES))
        for f in FAILURES[:12]:
            print("  " + f)
        return 1
    print("PASS: %d derived supplementary sheets reproduce from their records"
          % len(CHECKED))
    for c in CHECKED:
        print("  " + c)
    print("  OUT OF SCOPE (curated, not derived): "
          "whole_variant_negatives, benchmark_variants, model_scope, index")
    return 0


if __name__ == "__main__":
    sys.exit(main())
