#!/usr/bin/env python3
"""Verify canonical-to-ledger synchronization and metadata completeness.

This is a deliberately narrow synchronization/completeness gate. It proves
that each frozen canonical commitment has exactly one matching evidence-ledger
row, with a matching expected token and complete controlled metadata. It does
not independently adjudicate the literature or the scientific validity of any
label; that remains a separately governed source audit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Iterable


REPO = Path(__file__).resolve().parent.parent
DEFAULT_CANONICAL = REPO / "reference_outputs" / "scored_61var_canonical.csv"
DEFAULT_LEDGER = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"
DEFAULT_SUMMARY = REPO / "reference_outputs" / "COMAVI_evidence_ledger_summary.json"
# Re-pinned at v7.8 (BRCT provenance restoration).
#
# What authorised the move: the 12 brca1_brct rows carried NULL `gene` and NULL
# `structure_source` -- 20% of the benchmark and the entire monomer-fold arm,
# and the only rows in the table missing either column. Both values were
# already in the repository (inputs/raw/benchmark_variants_v6.csv: brca1,
# 1JNX.pdb) and had been dropped between the raw input and the canonical.
# scripts/apply_brct_provenance_v78.py fills ONLY null cells the raw input
# supplies, never overwrites a populated cell, and aborts if any metric column
# moves. No metric reads either column, and the headline numbers are
# bit-identical after the fill: mechanism consistency 0.6930 (39.5/57),
# all-row structural agreement 89/125.
#
# Nothing detected this gap: no gate reads these two columns. It surfaced only
# when building the per-variant table reviewers asked for.
#
# Previous pin (v7.7, axis_signature propagation):
#   85ed9f33d182d5c3c338106968d8076ab79883b21f3b930c935dc183271e8c7c
#
# --- v7.7 rationale, retained ---
#
# What authorised the move: the v7.6 ledger correction updated
# expected_mech_class on six variants but never propagated to axis_signature,
# the derived direction-explicit restatement of the same information that
# figures/src/figure3_axis_competency.py keys its class panel on. That panel
# therefore plotted a partition its own canonical contradicted.
# scripts/apply_axis_signature_v77.py repairs only rows where axis_signature
# CONTRADICTS expected_mech_class (exactly those six; --check passes on all 61),
# leaving the fourteen pre-existing drift rows from v7.2-v7.6 untouched. Headline
# metrics are unchanged by the repair: mechanism consistency 0.6930, all-row
# structural agreement 89/125. The invariant is now enforced by assertion in the
# consuming figure rather than asserted in a comment.
#
#
# v7.9 AUTHORISED RE-PIN. scripts/apply_ledger_corrections_v79.py reversed three
# binding tokens (msh2_msh6 C697F, kras_craf G12D, kras_craf G12V) from destab
# to neutral. Those three were the only rows in the ledger whose expected_token
# contradicted its own recorded evidence_basis, which all three state as intact
# binding; they were introduced by one v7.6 review batch (units D09, D18) that
# changed the tokens without changing the bases. The contradiction is now
# enforced by scripts/verify_ledger_token_basis.py rather than left to audit.
# Headline metrics DO move and the manuscript moves with them: mechanism
# consistency 0.6930 -> 0.7193, primary structural agreement 89/124 -> 92/124.
# This pin is updated deliberately, not to silence a failure -- the mismatch it
# raised was correct and is the only automatic guard on the canonical changing.
#
# Previous pins:
#   v7.8  fb80b9849d8bdee3ae9b511abf2e6da9656071208eacff5be08e23366ce2d614
#   v7.6  e4d657dee2625580bd41da05d9251fb69a224ddefdff6286f82a236511465d28
#   v7.5  88cee917d00ea6705e851b59b7551ef8211052011768a732462ee59ef45031bb
CANONICAL_SHA256 = "aba940a0ad211535b9d6e1d6ce83ca480b0dd4f38b1e14a033d183a9e1833e12"

# Committed-axis count. v7.6 withdrew 11 commitments (expected token -> unknown),
# taking the ledger and the canonical from 109 committed axes to 98.
COMMITTED_AXIS_COUNT = 98

AXES = {
    "monomer": "expected_ddg_monomer",
    "fold_complex": "expected_ddg_fold_complex",
    "binding": "expected_ddg_binding",
}
COMMITTED = {"destab", "neutral", "stab"}
EVIDENCE_TYPES = {
    "E1_quantitative_energetic",
    "E2_quantitative_functional",
    "E3_qualitative_experimental",
    "E4_population_frequency",
    "E5_inferred_no_axis_assay",
}
DIRECTNESS = {"direct", "coupled", "off_axis", "inferred"}
REQUIRED_LEDGER_COLUMNS = {
    "system",
    "variant",
    "axis",
    "expected_token",
    "evidence_type",
    "evidence_directness",
    "evidence_basis",
    "evidence_citation",
}
EXPECTED_A636P = {
    ("msh2_msh6", "A636P", "monomer"): {
        "expected_token": "neutral",
        "evidence_type": "E3_qualitative_experimental",
        "evidence_directness": "direct",
        "evidence_basis": (
            "MSH2 expression/stability preserved; ATPase mismatch "
            "binding/release defect"
        ),
        "evidence_citation": (
            "Ollila 2008 Hum Mutat 29:1355 (PMID 18951462); "
            "Ollila 2006 Gastroenterology (PMID 17101317)"
        ),
    },
    ("msh2_msh6", "A636P", "fold_complex"): {
        "expected_token": "neutral",
        "evidence_type": "E3_qualitative_experimental",
        "evidence_directness": "direct",
        "evidence_basis": (
            "MSH2-MSH6 interaction intact; ATPase mismatch binding/release defect"
        ),
        "evidence_citation": (
            "Ollila 2008 Hum Mutat 29:1355 (PMID 18951462); "
            "Ollila 2006 Gastroenterology (PMID 17101317)"
        ),
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def canonical_commitments(
    rows: Iterable[dict[str, str]],
) -> dict[tuple[str, str, str], str]:
    commitments: dict[tuple[str, str, str], str] = {}
    seen_variants: set[tuple[str, str]] = set()
    for row in rows:
        variant_key = (row["system"].strip(), row["variant"].strip())
        if variant_key in seen_variants:
            raise ValueError(f"duplicate canonical variant key: {variant_key}")
        seen_variants.add(variant_key)
        for axis, column in AXES.items():
            token = row.get(column, "").strip().lower()
            if token in COMMITTED:
                commitments[(variant_key[0], variant_key[1], axis)] = token
    return commitments


def ledger_map(
    rows: Iterable[dict[str, str]],
) -> dict[tuple[str, str, str], dict[str, str]]:
    ledger: dict[tuple[str, str, str], dict[str, str]] = {}
    for source_row in rows:
        row = {key: (value or "").strip() for key, value in source_row.items()}
        key = (row["system"], row["variant"], row["axis"])
        if key in ledger:
            raise ValueError(f"duplicate ledger key: {key}")
        ledger[key] = row
    return ledger


def expected_summary(
    ledger: dict[tuple[str, str, str], dict[str, str]],
) -> dict[str, object]:
    by_type = Counter(row["evidence_type"] for row in ledger.values())
    by_directness = Counter(
        row["evidence_directness"] for row in ledger.values()
    )
    by_axis = Counter(key[2] for key in ledger)
    type_by_axis = {
        axis: dict(
            Counter(
                row["evidence_type"]
                for key, row in ledger.items()
                if key[2] == axis
            )
        )
        for axis in AXES
    }
    count = len(ledger)
    energetic = by_type["E1_quantitative_energetic"]
    quantitative = energetic + by_type["E2_quantitative_functional"]
    return {
        "n_committed_axes": count,
        "n_variants": len({(key[0], key[1]) for key in ledger}),
        "n_systems": len({key[0] for key in ledger}),
        "by_evidence_type": dict(by_type),
        "by_directness": dict(by_directness),
        "by_axis": dict(by_axis),
        "n_energetic_direct": sum(
            1
            for row in ledger.values()
            if row["evidence_type"] == "E1_quantitative_energetic"
            and row["evidence_directness"] == "direct"
        ),
        "frac_energetic": round(energetic / count, 4),
        "frac_quantitative": round(quantitative / count, 4),
        "frac_direct": round(by_directness["direct"] / count, 4),
        "frac_inferred": round(
            by_type["E5_inferred_no_axis_assay"] / count, 4
        ),
        "type_by_axis": type_by_axis,
    }


def verify(canonical: Path, ledger: Path, summary: Path) -> list[str]:
    failures: list[str] = []

    observed_hash = sha256(canonical)
    if observed_hash != CANONICAL_SHA256:
        failures.append(
            f"canonical SHA-256 mismatch: {observed_hash} != {CANONICAL_SHA256}"
        )

    canonical_fields, canonical_rows = read_csv(canonical)
    missing_canonical = {"system", "variant", *AXES.values()} - set(
        canonical_fields
    )
    if missing_canonical:
        failures.append(f"canonical missing columns: {sorted(missing_canonical)}")
        return failures

    ledger_fields, ledger_rows = read_csv(ledger)
    missing_ledger = REQUIRED_LEDGER_COLUMNS - set(ledger_fields)
    if missing_ledger:
        failures.append(f"ledger missing columns: {sorted(missing_ledger)}")
        return failures

    try:
        canonical_map = canonical_commitments(canonical_rows)
    except ValueError as error:
        failures.append(str(error))
        canonical_map = None

    try:
        evidence_map = ledger_map(ledger_rows)
    except ValueError as error:
        failures.append(str(error))
        evidence_map = None

    if canonical_map is not None:
        if len(canonical_map) != COMMITTED_AXIS_COUNT:
            failures.append(
                f"canonical committed-axis count {len(canonical_map)} "
                f"!= {COMMITTED_AXIS_COUNT}"
            )

    if evidence_map is not None:
        if len(evidence_map) != COMMITTED_AXIS_COUNT:
            failures.append(
                f"ledger record count {len(evidence_map)} != {COMMITTED_AXIS_COUNT}"
            )

        for key, row in sorted(evidence_map.items()):
            if row["expected_token"] not in COMMITTED:
                failures.append(
                    f"invalid expected_token for {key}: {row['expected_token']!r}"
                )
            if row["evidence_type"] not in EVIDENCE_TYPES:
                failures.append(
                    f"invalid evidence_type for {key}: {row['evidence_type']!r}"
                )
            if row["evidence_directness"] not in DIRECTNESS:
                failures.append(
                    f"invalid evidence_directness for {key}: "
                    f"{row['evidence_directness']!r}"
                )
            for column in ("evidence_basis", "evidence_citation"):
                if not row[column]:
                    failures.append(f"blank {column} for {key}")

        for key, expected in EXPECTED_A636P.items():
            row = evidence_map.get(key)
            if row is None:
                failures.append(f"missing frozen A636P reconciliation row: {key}")
                continue
            for field, expected_value in expected.items():
                if row.get(field) != expected_value:
                    failures.append(
                        f"A636P {key[2]} {field} mismatch: "
                        f"{row.get(field)!r} != {expected_value!r}"
                    )

        if summary.exists():
            try:
                observed_summary = json.loads(summary.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                failures.append(f"summary is not valid JSON: {error}")
            else:
                if not isinstance(observed_summary, dict):
                    failures.append("summary root must be a JSON object")
                else:
                    for key, expected_value in expected_summary(
                        evidence_map
                    ).items():
                        if observed_summary.get(key) != expected_value:
                            failures.append(
                                f"summary field {key!r} mismatch: "
                                f"observed={observed_summary.get(key)!r}, "
                                f"expected={expected_value!r}"
                            )
        else:
            failures.append(f"summary file missing: {summary}")

    if canonical_map is not None and evidence_map is not None:
        missing = sorted(set(canonical_map) - set(evidence_map))
        extra = sorted(set(evidence_map) - set(canonical_map))
        if missing:
            failures.append(f"ledger missing {len(missing)} key(s): {missing}")
        if extra:
            failures.append(f"ledger has {len(extra)} extra key(s): {extra}")
        for key in sorted(set(canonical_map) & set(evidence_map)):
            observed_token = evidence_map[key]["expected_token"]
            if observed_token != canonical_map[key]:
                failures.append(
                    f"expected-token mismatch for {key}: "
                    f"ledger={observed_token!r}, canonical={canonical_map[key]!r}"
                )

    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, default=DEFAULT_CANONICAL)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    args = parser.parse_args()

    failures = verify(args.canonical, args.ledger, args.summary)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)

    _, ledger_rows = read_csv(args.ledger)
    counts = Counter(row["axis"] for row in ledger_rows)
    print(
        "PASS evidence ledger synchronization/completeness: "
        f"{len(ledger_rows)}/{COMMITTED_AXIS_COUNT} canonical commitments; "
        f"monomer={counts['monomer']}, "
        f"fold_complex={counts['fold_complex']}, "
        f"binding={counts['binding']}; "
        f"canonical_sha256={CANONICAL_SHA256}; "
        "scientific_adjudication=no"
    )


if __name__ == "__main__":
    main()
