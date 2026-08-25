#!/usr/bin/env python3
"""Sequence-aware setup helpers for the public COMAVI notebook.

The module deliberately separates *safe automation* from biological judgment.
It can align reference sequences to coordinate chains, infer uniform numbering
offsets, rank chain assignments, summarize variant coverage, and package a
validated system for reuse.  It does not decide whether a structure represents
the biologically relevant ligand, oligomer, conformation, or disease state.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import hashlib
import html
import itertools
import json
import math
import re
import shutil
import zipfile

import gemmi
import numpy as np
import pandas as pd
import yaml
from Bio import Align
from scipy.optimize import linear_sum_assignment

SETUP_WIZARD_VERSION = "COMAVI-Setup-v1"
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
STATUS_ORDER = {"green": 0, "yellow": 1, "red": 2}
STATUS_SYMBOL = {"green": "GREEN", "yellow": "YELLOW", "red": "RED"}


@dataclass(frozen=True)
class ChainRecord:
    chain_id: str
    sequence: str
    residue_numbers: tuple[int, ...]
    insertion_codes: tuple[str, ...]
    residue_names: tuple[str, ...]
    median_bfactor: float | None


@dataclass(frozen=True)
class AlignmentResult:
    gene: str
    chain_id: str
    identity: float
    reference_coverage: float
    chain_coverage: float
    aligned_pairs: int
    matches: int
    score: float
    reference_to_structure: dict[int, tuple[int, str, str]]
    numbering_offset: int | None
    numbering_uniformity: float


@dataclass(frozen=True)
class OffsetInference:
    gene: str
    selected_offset: int | None
    candidate_offsets: tuple[int, ...]
    matched_variants: int
    total_variants: int
    status: str
    message: str


@dataclass(frozen=True)
class GeneChainAssignment:
    gene: str
    chain_id: str | None
    identity: float | None
    reference_coverage: float | None
    chain_coverage: float | None
    variant_coverage: float
    assignment_margin: float | None
    raw_numbering_offset: int | None
    pipeline_offset: int | None
    numbering_uniformity: float | None
    status: str
    message: str


@dataclass(frozen=True)
class StructureAssessment:
    structure_path: str
    structure_name: str
    source_label: str
    status: str
    score: float
    chain_map: dict[str, str]
    pipeline_offsets: dict[str, int]
    input_offsets: dict[str, int]
    assignment_rows: tuple[GeneChainAssignment, ...]
    variant_rows: tuple[dict[str, Any], ...]
    messages: tuple[str, ...]


@dataclass(frozen=True)
class MonomerAssessment:
    gene: str
    structure_path: str
    chain_id: str | None
    input_offset: int | None
    pipeline_offset: int | None
    identity: float | None
    reference_coverage: float | None
    variant_coverage: float
    numbering_uniformity: float | None
    status: str
    message: str


def normalize_gene(value: Any) -> str:
    return str(value).strip().lower()


def normalize_sequence(value: Any) -> str:
    sequence = re.sub(r"[^A-Za-z]", "", str(value or "")).upper()
    return sequence


def _status_max(*statuses: str) -> str:
    values = [status for status in statuses if status in STATUS_ORDER]
    if not values:
        return "red"
    return max(values, key=lambda value: STATUS_ORDER[value])


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a simple FASTA file into normalized, lower-case gene keys."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []

    def commit() -> None:
        nonlocal current, chunks
        if current is None:
            return
        sequence = normalize_sequence("".join(chunks))
        if not sequence:
            raise ValueError(f"FASTA record {current!r} has no sequence in {path}.")
        if current in records and records[current] != sequence:
            raise ValueError(f"Conflicting FASTA records for {current!r} in {path}.")
        records[current] = sequence

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            commit()
            label = line[1:].strip().split()[0] if line[1:].strip() else ""
            if not label:
                raise ValueError(f"FASTA header is empty in {path}.")
            current = normalize_gene(label)
            chunks = []
        else:
            if current is None:
                raise ValueError(f"FASTA sequence appears before its header in {path}.")
            chunks.append(line)
    commit()
    if not records:
        raise ValueError(f"No FASTA records were found in {path}.")
    return records


def write_fasta(path: str | Path, records: Mapping[str, str]) -> Path:
    """Write normalized reference sequences deterministically."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for gene, sequence in sorted(records.items()):
        normalized = normalize_sequence(sequence)
        if not normalized:
            raise ValueError(f"Reference sequence is empty for {gene!r}.")
        lines.extend([f">{normalize_gene(gene)}", normalized])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _variant_frame(variants: pd.DataFrame) -> pd.DataFrame:
    required = ["gene", "ref_aa", "position", "alt_aa"]
    missing = [column for column in required if column not in variants.columns]
    if missing:
        raise ValueError("Variant table lacks required columns: " + ", ".join(missing))
    frame = variants.copy()
    frame["gene"] = frame["gene"].map(normalize_gene)
    frame["ref_aa"] = frame["ref_aa"].astype(str).str.strip().str.upper()
    frame["alt_aa"] = frame["alt_aa"].astype(str).str.strip().str.upper()
    frame["position"] = pd.to_numeric(frame["position"], errors="raise").astype(int)
    if "variant" not in frame.columns:
        frame["variant"] = (
            frame["ref_aa"] + frame["position"].astype(str) + frame["alt_aa"]
        )
    return frame


def parse_gene_integer_map(
    raw_text: Any,
    *,
    configured_genes: Iterable[str] | None = None,
    label: str = "override",
) -> dict[str, int]:
    """Parse comma, semicolon, or newline separated ``GENE:INTEGER`` pairs."""
    text = str(raw_text or "").strip()
    if not text:
        return {}
    configured = (
        {normalize_gene(value) for value in configured_genes}
        if configured_genes is not None
        else None
    )
    result: dict[str, int] = {}
    for token in re.split(r"[,;\n]+", text):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            raise ValueError(f"{label} entry {token!r} must use GENE:INTEGER.")
        gene, raw_value = (part.strip() for part in token.split(":", 1))
        gene_key = normalize_gene(gene)
        if configured is not None and gene_key not in configured:
            raise ValueError(
                f"{label} specifies {gene!r}, which is not in the configured genes: "
                + ", ".join(sorted(configured))
            )
        try:
            value = int(raw_value)
        except ValueError as error:
            raise ValueError(
                f"{label} value for {gene!r} must be an integer; received {raw_value!r}."
            ) from error
        if gene_key in result and result[gene_key] != value:
            raise ValueError(
                f"{label} gives conflicting values for {gene!r}: "
                f"{result[gene_key]} and {value}."
            )
        result[gene_key] = value
    return result


def parse_gene_text_map(
    raw_text: Any,
    *,
    configured_genes: Iterable[str] | None = None,
    label: str = "chain map",
) -> dict[str, str]:
    """Parse comma, semicolon, or newline separated ``GENE:TEXT`` pairs."""
    text = str(raw_text or "").strip()
    if not text:
        return {}
    configured = (
        {normalize_gene(value) for value in configured_genes}
        if configured_genes is not None
        else None
    )
    result: dict[str, str] = {}
    for token in re.split(r"[,;\n]+", text):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            raise ValueError(f"{label} entry {token!r} must use GENE:CHAIN.")
        gene, raw_value = (part.strip() for part in token.split(":", 1))
        gene_key = normalize_gene(gene)
        if configured is not None and gene_key not in configured:
            raise ValueError(
                f"{label} specifies {gene!r}, which is not in the configured genes: "
                + ", ".join(sorted(configured))
            )
        if not raw_value:
            raise ValueError(f"{label} gives an empty chain for {gene!r}.")
        if gene_key in result and result[gene_key] != raw_value:
            raise ValueError(
                f"{label} gives conflicting chains for {gene!r}: "
                f"{result[gene_key]!r} and {raw_value!r}."
            )
        result[gene_key] = raw_value
    return result


def infer_input_offset(
    gene: str,
    variants: pd.DataFrame,
    reference_sequence: str,
    *,
    override: int | None = None,
) -> OffsetInference:
    """Infer the uniform submitted-numbering offset relative to a reference sequence.

    COMAVI uses ``structure_position = submitted_position - offset``.  Here the
    same convention is applied to the reference sequence: ``reference_position
    = submitted_position - selected_offset``.
    """
    gene_key = normalize_gene(gene)
    frame = _variant_frame(variants)
    frame = frame.loc[frame["gene"].eq(gene_key)].copy()
    sequence = normalize_sequence(reference_sequence)
    if not sequence:
        return OffsetInference(
            gene_key,
            None,
            (),
            0,
            len(frame),
            "red",
            "No reference sequence is available.",
        )
    if frame.empty:
        return OffsetInference(
            gene_key,
            int(override or 0),
            (int(override or 0),),
            0,
            0,
            "green",
            "No submitted variants belong to this gene; canonical input numbering (offset 0) is assumed and future variants will be rechecked.",
        )

    candidate_sets: list[set[int]] = []
    for row in frame.itertuples(index=False):
        positions = {
            index + 1
            for index, amino_acid in enumerate(sequence)
            if amino_acid == str(row.ref_aa).upper()
        }
        candidate_sets.append({int(row.position) - position for position in positions})

    intersection = set.intersection(*candidate_sets) if candidate_sets else set()
    all_candidates = set.union(*candidate_sets) if candidate_sets else set()

    def match_count(offset: int) -> int:
        total = 0
        for row in frame.itertuples(index=False):
            reference_position = int(row.position) - offset
            if (
                1 <= reference_position <= len(sequence)
                and sequence[reference_position - 1] == str(row.ref_aa).upper()
            ):
                total += 1
        return total

    if override is not None:
        selected = int(override)
        matches = match_count(selected)
        status = "green" if matches == len(frame) else "red"
        message = (
            f"Manual input-numbering offset {selected} matches {matches}/{len(frame)} variants."
        )
        return OffsetInference(
            gene_key,
            selected,
            tuple(sorted(intersection, key=lambda value: (abs(value), value))),
            matches,
            len(frame),
            status,
            message,
        )

    if intersection:
        ordered = tuple(sorted(intersection, key=lambda value: (abs(value), value)))
        selected = ordered[0]
        status = "green" if len(ordered) == 1 else "yellow"
        if len(ordered) == 1:
            message = f"A unique numbering offset ({selected}) matches all submitted variants."
        else:
            message = (
                f"Multiple offsets match all variants; selected the smallest absolute offset "
                f"({selected}). Alternatives: {list(ordered[1:8])}. Confirm before reuse."
            )
        return OffsetInference(
            gene_key,
            selected,
            ordered,
            len(frame),
            len(frame),
            status,
            message,
        )

    if not all_candidates:
        return OffsetInference(
            gene_key,
            None,
            (),
            0,
            len(frame),
            "red",
            "No reference-sequence position carries the submitted reference amino acid.",
        )

    ranked = sorted(
        all_candidates,
        key=lambda value: (-match_count(value), abs(value), value),
    )
    selected = ranked[0]
    matches = match_count(selected)
    return OffsetInference(
        gene_key,
        selected,
        tuple(ranked[:20]),
        matches,
        len(frame),
        "red",
        (
            f"No single offset matches every submitted variant. Best candidate {selected} "
            f"matches {matches}/{len(frame)}; check isoform and variant nomenclature."
        ),
    )


def _aa_from_residue_name(name: str) -> str | None:
    info = gemmi.find_tabulated_residue(str(name).strip().upper())
    if info.kind != gemmi.ResidueKind.AA:
        return None
    code = str(info.one_letter_code or "X")
    if code == "m":
        code = "M"
    code = code.upper()
    return code if len(code) == 1 else "X"


def extract_chain_records(structure_path: str | Path) -> dict[str, ChainRecord]:
    """Extract observed amino-acid sequences and residue identifiers by chain."""
    path = Path(structure_path)
    if not path.is_file():
        raise FileNotFoundError(path)
    structure = gemmi.read_structure(str(path))
    if len(structure) == 0:
        raise ValueError(f"Structure contains no models: {path}")
    model = structure[0]
    records: dict[str, ChainRecord] = {}
    for chain in model:
        sequence: list[str] = []
        residue_numbers: list[int] = []
        insertion_codes: list[str] = []
        residue_names: list[str] = []
        b_factors: list[float] = []
        for residue in chain:
            amino_acid = _aa_from_residue_name(residue.name)
            if amino_acid is None:
                continue
            atoms = list(residue)
            if not atoms:
                continue
            sequence.append(amino_acid)
            residue_numbers.append(int(residue.seqid.num))
            insertion_codes.append(str(residue.seqid.icode or "").strip())
            residue_names.append(str(residue.name).strip())
            finite_b = [float(atom.b_iso) for atom in atoms if math.isfinite(float(atom.b_iso))]
            if finite_b:
                b_factors.append(float(np.median(finite_b)))
        if sequence:
            chain_id = str(chain.name)
            records[chain_id] = ChainRecord(
                chain_id=chain_id,
                sequence="".join(sequence),
                residue_numbers=tuple(residue_numbers),
                insertion_codes=tuple(insertion_codes),
                residue_names=tuple(residue_names),
                median_bfactor=(float(np.median(b_factors)) if b_factors else None),
            )
    if not records:
        raise ValueError(f"No amino-acid chains were found in {path}")
    return records


def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5
    # Free terminal gaps are appropriate for domain constructs and unresolved ends.
    try:
        aligner.end_insertion_score = 0.0
        aligner.end_deletion_score = 0.0
    except AttributeError:
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
    return aligner


def align_reference_to_chain(
    gene: str,
    reference_sequence: str,
    chain: ChainRecord,
) -> AlignmentResult:
    reference = normalize_sequence(reference_sequence)
    if not reference:
        raise ValueError(f"Reference sequence is empty for {gene}")
    if not chain.sequence:
        raise ValueError(f"Chain {chain.chain_id} has no amino-acid sequence")
    alignment = _aligner().align(reference, chain.sequence)[0]
    mapping: dict[int, tuple[int, str, str]] = {}
    matches = 0
    aligned_pairs = 0
    reference_segments, chain_segments = alignment.aligned
    for (reference_start, reference_end), (chain_start, chain_end) in zip(
        reference_segments, chain_segments
    ):
        length = min(reference_end - reference_start, chain_end - chain_start)
        for index in range(length):
            reference_index = int(reference_start + index)
            chain_index = int(chain_start + index)
            reference_position = reference_index + 1
            mapping[reference_position] = (
                int(chain.residue_numbers[chain_index]),
                chain.insertion_codes[chain_index],
                chain.sequence[chain_index],
            )
            aligned_pairs += 1
            if reference[reference_index] == chain.sequence[chain_index]:
                matches += 1
    identity = matches / aligned_pairs if aligned_pairs else 0.0
    reference_coverage = aligned_pairs / len(reference) if reference else 0.0
    chain_coverage = aligned_pairs / len(chain.sequence) if chain.sequence else 0.0

    raw_offsets = [
        reference_position - structure_info[0]
        for reference_position, structure_info in mapping.items()
        if not structure_info[1]
    ]
    numbering_offset = None
    numbering_uniformity = 0.0
    if raw_offsets:
        values, counts = np.unique(np.asarray(raw_offsets, dtype=int), return_counts=True)
        best_index = int(np.argmax(counts))
        numbering_offset = int(values[best_index])
        numbering_uniformity = float(counts[best_index] / len(raw_offsets))

    return AlignmentResult(
        gene=normalize_gene(gene),
        chain_id=chain.chain_id,
        identity=float(identity),
        reference_coverage=float(reference_coverage),
        chain_coverage=float(chain_coverage),
        aligned_pairs=int(aligned_pairs),
        matches=int(matches),
        score=float(alignment.score),
        reference_to_structure=mapping,
        numbering_offset=numbering_offset,
        numbering_uniformity=numbering_uniformity,
    )


def _variant_reference_positions(
    gene: str,
    variants: pd.DataFrame,
    input_offset: int | None,
) -> list[tuple[int, str, str, int]]:
    frame = _variant_frame(variants)
    gene_key = normalize_gene(gene)
    rows = frame.loc[frame["gene"].eq(gene_key)]
    if input_offset is None:
        return []
    return [
        (
            int(row.position) - int(input_offset),
            str(row.ref_aa).upper(),
            str(row.variant),
            int(row.position),
        )
        for row in rows.itertuples(index=False)
    ]


def _alignment_variant_coverage(
    alignment: AlignmentResult,
    variant_positions: Sequence[tuple[int, str, str, int]],
) -> tuple[float, list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    covered = 0
    for reference_position, expected_aa, variant, submitted_position in variant_positions:
        mapped = alignment.reference_to_structure.get(reference_position)
        structure_position = mapped[0] if mapped else None
        insertion_code = mapped[1] if mapped else None
        structure_aa = mapped[2] if mapped else None
        passed = bool(
            mapped
            and not insertion_code
            and structure_aa == expected_aa
        )
        if passed:
            covered += 1
        rows.append(
            {
                "variant": variant,
                "submitted_position": submitted_position,
                "reference_position": reference_position,
                "expected_aa": expected_aa,
                "structure_position": structure_position,
                "insertion_code": insertion_code,
                "structure_aa": structure_aa,
                "covered": passed,
            }
        )
    coverage = covered / len(variant_positions) if variant_positions else 1.0
    return float(coverage), rows


def _alignment_rank_score(
    alignment: AlignmentResult,
    variant_coverage: float,
) -> float:
    return float(
        200.0 * variant_coverage
        + 100.0 * alignment.identity
        + 25.0 * min(alignment.reference_coverage, 1.0)
        + 10.0 * min(alignment.chain_coverage, 1.0)
        + 10.0 * alignment.numbering_uniformity
    )


def assess_structure(
    structure_path: str | Path,
    reference_sequences: Mapping[str, str],
    variants: pd.DataFrame,
    *,
    source_label: str = "structure",
    input_offset_overrides: Mapping[str, int] | None = None,
    chain_overrides: Mapping[str, str] | None = None,
    pipeline_offset_overrides: Mapping[str, int] | None = None,
) -> StructureAssessment:
    """Align configured genes to a structure and derive safe COMAVI mappings."""
    path = Path(structure_path)
    frame = _variant_frame(variants)
    references = {
        normalize_gene(gene): normalize_sequence(sequence)
        for gene, sequence in reference_sequences.items()
    }
    genes = list(references)
    if len(genes) != len(set(genes)):
        raise ValueError("Duplicate gene identifiers are not supported by this setup wizard version.")
    if not genes:
        raise ValueError("No reference sequences were supplied.")
    chains = extract_chain_records(path)
    input_offset_overrides = {
        normalize_gene(key): int(value)
        for key, value in (input_offset_overrides or {}).items()
    }
    chain_overrides = {
        normalize_gene(key): str(value)
        for key, value in (chain_overrides or {}).items()
    }
    pipeline_offset_overrides = {
        normalize_gene(key): int(value)
        for key, value in (pipeline_offset_overrides or {}).items()
    }

    offset_inferences: dict[str, OffsetInference] = {}
    alignments: dict[tuple[str, str], AlignmentResult] = {}
    variant_positions: dict[str, list[tuple[int, str, str, int]]] = {}
    variant_coverage: dict[tuple[str, str], float] = {}
    variant_detail: dict[tuple[str, str], list[dict[str, Any]]] = {}
    rank_score: dict[tuple[str, str], float] = {}

    for gene in genes:
        inference = infer_input_offset(
            gene,
            frame,
            references[gene],
            override=input_offset_overrides.get(gene),
        )
        offset_inferences[gene] = inference
        positions = _variant_reference_positions(gene, frame, inference.selected_offset)
        variant_positions[gene] = positions
        for chain_id, chain in chains.items():
            alignment = align_reference_to_chain(gene, references[gene], chain)
            coverage, details = _alignment_variant_coverage(alignment, positions)
            alignments[(gene, chain_id)] = alignment
            variant_coverage[(gene, chain_id)] = coverage
            variant_detail[(gene, chain_id)] = details
            rank_score[(gene, chain_id)] = _alignment_rank_score(alignment, coverage)

    assignments: dict[str, str] = {}
    used_chains: set[str] = set()
    for gene, chain_id in chain_overrides.items():
        if gene not in references:
            raise ValueError(f"Chain override refers to unconfigured gene {gene!r}.")
        if chain_id not in chains:
            raise ValueError(
                f"Chain override {gene}:{chain_id} is absent from {path.name}; "
                f"available chains: {sorted(chains)}"
            )
        if chain_id in used_chains:
            raise ValueError(f"Chain {chain_id!r} is assigned to more than one gene.")
        assignments[gene] = chain_id
        used_chains.add(chain_id)

    remaining_genes = [gene for gene in genes if gene not in assignments]
    remaining_chains = [chain_id for chain_id in chains if chain_id not in used_chains]
    if remaining_genes and remaining_chains:
        matrix = np.asarray(
            [
                [rank_score[(gene, chain_id)] for chain_id in remaining_chains]
                for gene in remaining_genes
            ],
            dtype=float,
        )
        rows, columns = linear_sum_assignment(-matrix)
        for row_index, column_index in zip(rows, columns):
            assignments[remaining_genes[int(row_index)]] = remaining_chains[int(column_index)]

    assignment_rows: list[GeneChainAssignment] = []
    all_variant_rows: list[dict[str, Any]] = []
    chain_map: dict[str, str] = {}
    pipeline_offsets: dict[str, int] = {}
    messages: list[str] = []

    for gene in genes:
        inference = offset_inferences[gene]
        chain_id = assignments.get(gene)
        if chain_id is None:
            row = GeneChainAssignment(
                gene=gene,
                chain_id=None,
                identity=None,
                reference_coverage=None,
                chain_coverage=None,
                variant_coverage=0.0,
                assignment_margin=None,
                raw_numbering_offset=None,
                pipeline_offset=None,
                numbering_uniformity=None,
                status="red",
                message="No unassigned amino-acid chain was available for this gene.",
            )
            assignment_rows.append(row)
            messages.append(f"{gene}: {row.message}")
            continue

        alignment = alignments[(gene, chain_id)]
        coverage = variant_coverage[(gene, chain_id)]
        alternatives = sorted(
            (
                rank_score[(gene, other_chain)]
                for other_chain in chains
                if other_chain != chain_id
            ),
            reverse=True,
        )
        margin = (
            rank_score[(gene, chain_id)] - alternatives[0]
            if alternatives
            else None
        )
        raw_offset = alignment.numbering_offset
        pipeline_offset = None
        if inference.selected_offset is not None and raw_offset is not None:
            pipeline_offset = int(inference.selected_offset + raw_offset)
        if gene in pipeline_offset_overrides:
            pipeline_offset = pipeline_offset_overrides[gene]

        status = "green"
        reasons: list[str] = []
        if inference.status == "red":
            status = "red"
            reasons.append(inference.message)
        elif inference.status == "yellow":
            status = "yellow"
            reasons.append(inference.message)
        if alignment.identity < 0.70:
            status = "red"
            reasons.append(f"chain identity is only {alignment.identity:.1%}")
        elif alignment.identity < 0.90 and status != "red":
            status = "yellow"
            reasons.append(f"chain identity is {alignment.identity:.1%}")
        if coverage < 1.0:
            status = "red"
            reasons.append(f"only {coverage:.0%} of submitted variants map with the expected residue")
        if alignment.numbering_uniformity < 0.95:
            status = "red"
            reasons.append(
                "structure numbering is not represented by one uniform offset; use a preprocessed structure"
            )
        elif alignment.numbering_uniformity < 1.0 and status != "red":
            status = "yellow"
            reasons.append(
                f"numbering is {alignment.numbering_uniformity:.1%} uniform across aligned residues"
            )
        if margin is not None and margin < 15.0 and status == "green":
            status = "yellow"
            reasons.append("another chain has a similar sequence-assignment score")
        if pipeline_offset is None:
            status = "red"
            reasons.append("a pipeline numbering offset could not be derived")

        if not reasons:
            reasons.append(
                f"chain {chain_id} maps at {alignment.identity:.1%} identity and covers every submitted variant"
            )
        message = "; ".join(reasons)
        assignment_rows.append(
            GeneChainAssignment(
                gene=gene,
                chain_id=chain_id,
                identity=alignment.identity,
                reference_coverage=alignment.reference_coverage,
                chain_coverage=alignment.chain_coverage,
                variant_coverage=coverage,
                assignment_margin=margin,
                raw_numbering_offset=raw_offset,
                pipeline_offset=pipeline_offset,
                numbering_uniformity=alignment.numbering_uniformity,
                status=status,
                message=message,
            )
        )
        chain_map[gene] = chain_id
        if pipeline_offset is not None:
            pipeline_offsets[gene] = int(pipeline_offset)
        messages.append(f"{gene}: {message}")

        for detail in variant_detail[(gene, chain_id)]:
            detail = dict(detail)
            detail.update(
                {
                    "gene": gene,
                    "chain_id": chain_id,
                    "input_offset": inference.selected_offset,
                    "pipeline_offset": pipeline_offset,
                    "status": "green" if detail["covered"] else "red",
                }
            )
            all_variant_rows.append(detail)

    overall_status = "green"
    for row in assignment_rows:
        overall_status = _status_max(overall_status, row.status)
    score = 0.0
    if assignment_rows:
        green_fraction = sum(row.status == "green" for row in assignment_rows) / len(assignment_rows)
        yellow_fraction = sum(row.status == "yellow" for row in assignment_rows) / len(assignment_rows)
        mean_identity = np.mean([row.identity or 0.0 for row in assignment_rows])
        mean_variant_coverage = np.mean([row.variant_coverage for row in assignment_rows])
        score = float(
            500.0 * green_fraction
            + 250.0 * yellow_fraction
            + 100.0 * mean_identity
            + 200.0 * mean_variant_coverage
        )

    return StructureAssessment(
        structure_path=str(path),
        structure_name=path.name,
        source_label=str(source_label),
        status=overall_status,
        score=score,
        chain_map=chain_map,
        pipeline_offsets=pipeline_offsets,
        input_offsets={
            gene: int(inference.selected_offset)
            for gene, inference in offset_inferences.items()
            if inference.selected_offset is not None
        },
        assignment_rows=tuple(assignment_rows),
        variant_rows=tuple(all_variant_rows),
        messages=tuple(messages),
    )


def assess_monomer(
    gene: str,
    structure_path: str | Path,
    reference_sequence: str,
    variants: pd.DataFrame,
    *,
    input_offset_override: int | None = None,
    chain_override: str | None = None,
    pipeline_offset_override: int | None = None,
) -> MonomerAssessment:
    gene_key = normalize_gene(gene)
    assessment = assess_structure(
        structure_path,
        {gene_key: reference_sequence},
        variants,
        source_label="monomer",
        input_offset_overrides=(
            {gene_key: input_offset_override}
            if input_offset_override is not None
            else None
        ),
        chain_overrides=({gene_key: chain_override} if chain_override else None),
        pipeline_offset_overrides=(
            {gene_key: pipeline_offset_override}
            if pipeline_offset_override is not None
            else None
        ),
    )
    row = assessment.assignment_rows[0]
    return MonomerAssessment(
        gene=gene_key,
        structure_path=str(structure_path),
        chain_id=row.chain_id,
        input_offset=assessment.input_offsets.get(gene_key),
        pipeline_offset=row.pipeline_offset,
        identity=row.identity,
        reference_coverage=row.reference_coverage,
        variant_coverage=row.variant_coverage,
        numbering_uniformity=row.numbering_uniformity,
        status=row.status,
        message=row.message,
    )


def assess_monomer_set(
    monomer_paths: Mapping[str, str | Path],
    reference_sequences: Mapping[str, str],
    variants: pd.DataFrame,
    *,
    input_offset_overrides: Mapping[str, int] | None = None,
    pipeline_offset_overrides: Mapping[str, int] | None = None,
) -> dict[str, MonomerAssessment]:
    results: dict[str, MonomerAssessment] = {}
    for gene, path in monomer_paths.items():
        gene_key = normalize_gene(gene)
        if gene_key not in reference_sequences:
            raise ValueError(f"Reference sequence is missing for monomer gene {gene_key}.")
        results[gene_key] = assess_monomer(
            gene_key,
            path,
            reference_sequences[gene_key],
            variants,
            input_offset_override=(input_offset_overrides or {}).get(gene_key),
            pipeline_offset_override=(pipeline_offset_overrides or {}).get(gene_key),
        )
    return results


def assessment_assignment_frame(assessment: StructureAssessment) -> pd.DataFrame:
    return pd.DataFrame([asdict(row) for row in assessment.assignment_rows])


def assessment_variant_frame(assessment: StructureAssessment) -> pd.DataFrame:
    return pd.DataFrame(list(assessment.variant_rows))


def monomer_assessment_frame(assessments: Mapping[str, MonomerAssessment]) -> pd.DataFrame:
    return pd.DataFrame([asdict(value) for value in assessments.values()])


def rank_structure_candidates(
    candidates: Sequence[Mapping[str, Any]],
    reference_sequences: Mapping[str, str],
    variants: pd.DataFrame,
    *,
    input_offset_overrides: Mapping[str, int] | None = None,
    maximum_candidates: int | None = None,
) -> tuple[pd.DataFrame, dict[str, StructureAssessment]]:
    """Assess local candidate coordinate files and return a ranked table.

    Each candidate mapping needs ``path`` and may include ``pdb_id``,
    ``resolution_A``, ``method``, ``title``, and ``year``.
    """
    rows: list[dict[str, Any]] = []
    assessments: dict[str, StructureAssessment] = {}
    iterable = candidates[:maximum_candidates] if maximum_candidates else candidates
    for index, candidate in enumerate(iterable, start=1):
        path = Path(candidate["path"])
        candidate_id = str(candidate.get("pdb_id") or path.stem).upper()
        try:
            assessment = assess_structure(
                path,
                reference_sequences,
                variants,
                source_label=candidate_id,
                input_offset_overrides=input_offset_overrides,
            )
            assessments[candidate_id] = assessment
            resolution = candidate.get("resolution_A")
            resolution_bonus = 0.0
            try:
                if resolution is not None and float(resolution) > 0:
                    resolution_bonus = max(0.0, 25.0 - 5.0 * float(resolution))
            except (TypeError, ValueError):
                pass
            rank_value = assessment.score + resolution_bonus
            rows.append(
                {
                    "rank_input": index,
                    "candidate_id": candidate_id,
                    "status": assessment.status,
                    "score": rank_value,
                    "chain_map": ", ".join(
                        f"{gene.upper()}:{chain}"
                        for gene, chain in sorted(assessment.chain_map.items())
                    ),
                    "pipeline_offsets": ", ".join(
                        f"{gene.upper()}:{offset}"
                        for gene, offset in sorted(assessment.pipeline_offsets.items())
                    ),
                    "resolution_A": resolution,
                    "method": candidate.get("method"),
                    "year": candidate.get("year"),
                    "title": candidate.get("title"),
                    "path": str(path),
                    "message": " | ".join(assessment.messages),
                }
            )
        except Exception as error:
            rows.append(
                {
                    "rank_input": index,
                    "candidate_id": candidate_id,
                    "status": "red",
                    "score": -1.0,
                    "chain_map": "",
                    "pipeline_offsets": "",
                    "resolution_A": candidate.get("resolution_A"),
                    "method": candidate.get("method"),
                    "year": candidate.get("year"),
                    "title": candidate.get("title"),
                    "path": str(path),
                    "message": f"Assessment failed: {error}",
                }
            )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        status_rank = frame["status"].map({"green": 0, "yellow": 1, "red": 2}).fillna(3)
        frame = (
            frame.assign(_status_rank=status_rank)
            .sort_values(["_status_rank", "score", "resolution_A"], ascending=[True, False, True], kind="mergesort", na_position="last")
            .drop(columns=["_status_rank"])
            .reset_index(drop=True)
        )
        frame.insert(0, "recommended_rank", np.arange(1, len(frame) + 1))
    return frame, assessments


def build_preflight_table(
    variants: pd.DataFrame,
    monomers: Mapping[str, MonomerAssessment],
    complex_assessment: StructureAssessment,
) -> pd.DataFrame:
    frame = _variant_frame(variants)
    complex_variants = assessment_variant_frame(complex_assessment)
    complex_lookup = {
        (normalize_gene(row["gene"]), str(row["variant"])): row
        for row in complex_variants.to_dict(orient="records")
    }
    rows: list[dict[str, Any]] = []
    for row in frame.itertuples(index=False):
        gene = normalize_gene(row.gene)
        variant = str(row.variant)
        monomer = monomers.get(gene)
        complex_row = complex_lookup.get((gene, variant))
        input_status = "green"
        monomer_status = monomer.status if monomer else "red"
        complex_status = (
            str(complex_row.get("status"))
            if complex_row is not None
            else "red"
        )
        overall = _status_max(input_status, monomer_status, complex_status)
        actions: list[str] = []
        if monomer is None:
            actions.append("Provide or remap a monomer structure.")
        elif monomer.status != "green":
            actions.append(monomer.message)
        if complex_row is None:
            actions.append("Map the variant to a complex chain.")
        elif complex_status != "green":
            actions.append("Check complex chain, numbering, and residue coverage.")
        if overall == "green":
            actions.append("Ready for COMAVI structural scoring.")
        rows.append(
            {
                "gene": gene,
                "variant": variant,
                "input": input_status,
                "monomer_mapping": monomer_status,
                "complex_mapping": complex_status,
                "overall": overall,
                "action": " ".join(actions),
                "monomer_chain": monomer.chain_id if monomer else None,
                "monomer_offset": monomer.pipeline_offset if monomer else None,
                "complex_chain": complex_row.get("chain_id") if complex_row else None,
                "complex_offset": complex_row.get("pipeline_offset") if complex_row else None,
                "complex_structure_position": complex_row.get("structure_position") if complex_row else None,
            }
        )
    return pd.DataFrame(rows)


def render_traffic_light_html(frame: pd.DataFrame, *, title: str = "COMAVI setup preflight") -> str:
    colors = {"green": "#d9ead3", "yellow": "#fff2cc", "red": "#f4cccc"}
    labels = {"green": "READY", "yellow": "REVIEW", "red": "STOP"}
    body: list[str] = []
    for row in frame.to_dict(orient="records"):
        status = str(row.get("overall", "red"))
        body.append(
            "<tr>"
            f"<td><b>{html.escape(str(row.get('gene', '')).upper())} "
            f"{html.escape(str(row.get('variant', '')))}</b></td>"
            f"<td style='background:{colors.get(status, colors['red'])};font-weight:bold'>"
            f"{labels.get(status, 'STOP')}</td>"
            f"<td>{html.escape(str(row.get('monomer_mapping', '')))}</td>"
            f"<td>{html.escape(str(row.get('complex_mapping', '')))}</td>"
            f"<td>{html.escape(str(row.get('action', '')))}</td>"
            "</tr>"
        )
    return (
        "<div style='font-family:Arial,sans-serif'>"
        f"<h3>{html.escape(title)}</h3>"
        "<table style='border-collapse:collapse;width:100%'>"
        "<thead><tr><th>Variant</th><th>Status</th><th>Monomer</th>"
        "<th>Complex</th><th>Next action</th></tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></div>"
    )


def _axis_explanation(axis: Any) -> tuple[str, str]:
    token = str(axis or "").strip().lower()
    if token in {"monomer", "monomer_fold", "fold_monomer"}:
        return (
            "isolated-protein stability",
            "Start with abundance, soluble-expression, thermal-shift, or chemical-stability measurements.",
        )
    if token in {"complex_context", "complex-context", "fold", "assembled_complex"}:
        return (
            "stability in the assembled complex",
            "Start with complex-assembly, co-expression, or subunit-stability measurements in the partner context.",
        )
    if token in {"binding", "interface"}:
        return (
            "partner interaction",
            "Start with affinity, pull-down, co-immunoprecipitation, or another direct interaction assay.",
        )
    return (
        "no single dominant modeled axis",
        "Inspect the complete signed mechanism profile and the model-scope limitations before choosing an assay.",
    )


def plain_language_summary_frame(results: pd.DataFrame) -> pd.DataFrame:
    frame = results.copy()
    if frame.empty:
        return pd.DataFrame()
    scores = pd.to_numeric(frame.get("isds_v1"), errors="coerce")
    ranks = scores.rank(method="min", ascending=False)
    total_available = int(scores.notna().sum())
    rows: list[dict[str, Any]] = []
    for index, row in frame.iterrows():
        axis_label, experiment = _axis_explanation(row.get("isds_dominant_axis"))
        score = scores.loc[index] if index in scores.index else np.nan
        rank = ranks.loc[index] if index in ranks.index else np.nan
        mechanism = None
        for column in ("comavi_mechanism_t25", "comavi_mechanism_t10", "comavi_mechanism"):
            if column in frame.columns and pd.notna(row.get(column)):
                mechanism = str(row.get(column))
                break
        rows.append(
            {
                "gene": str(row.get("gene", "")).upper(),
                "variant": str(row.get("variant", "")),
                "isds_v1": (float(score) if pd.notna(score) else None),
                "priority_rank": (int(rank) if pd.notna(rank) else None),
                "priority_denominator": total_available,
                "dominant_modeled_effect": axis_label,
                "mechanism_call": mechanism,
                "suggested_first_experiment": experiment,
                "interpretation_limit": (
                    "ISDS-v1 is a structural-prioritization index, not a probability or pathogenicity verdict."
                ),
            }
        )
    return pd.DataFrame(rows)


def render_plain_language_cards(results: pd.DataFrame) -> str:
    summary = plain_language_summary_frame(results)
    cards: list[str] = []
    for row in summary.to_dict(orient="records"):
        score = "unavailable" if row["isds_v1"] is None else f"{row['isds_v1']:.3f}"
        rank = (
            "not ranked"
            if row["priority_rank"] is None
            else f"{row['priority_rank']} of {row['priority_denominator']}"
        )
        cards.append(
            "<div style='border:1px solid #bbb;border-radius:8px;padding:12px;margin:8px 0'>"
            f"<h3 style='margin-top:0'>{html.escape(row['gene'])} {html.escape(row['variant'])}</h3>"
            f"<p><b>Structural-follow-up priority:</b> ISDS-v1 {score}; rank {rank}.</p>"
            f"<p><b>Dominant modeled effect:</b> {html.escape(row['dominant_modeled_effect'])}.</p>"
            f"<p><b>Mechanism call:</b> {html.escape(str(row['mechanism_call'] or 'no thresholded call shown'))}.</p>"
            f"<p><b>Suggested first experiment:</b> {html.escape(row['suggested_first_experiment'])}</p>"
            f"<p style='color:#555'><b>Important:</b> {html.escape(row['interpretation_limit'])}</p>"
            "</div>"
        )
    return "<div style='font-family:Arial,sans-serif'>" + "".join(cards) + "</div>"


def create_system_setup_bundle(
    destination: str | Path,
    *,
    config_path: str | Path,
    structure_paths: Sequence[str | Path],
    reference_sequences: Mapping[str, str],
    preflight: pd.DataFrame,
    setup_report: Mapping[str, Any],
) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    config_path = Path(config_path)
    structure_paths = [Path(path) for path in structure_paths]
    for path in [config_path, *structure_paths]:
        if not path.is_file():
            raise FileNotFoundError(path)
    manifest: dict[str, Any] = {
        "setup_wizard_version": SETUP_WIZARD_VERSION,
        "config": "config.yaml",
        "structures": [],
        "files": {},
    }
    temporary = destination.parent / (destination.stem + "_staging")
    if temporary.exists():
        shutil.rmtree(temporary)
    (temporary / "structures").mkdir(parents=True)
    shutil.copy2(config_path, temporary / "config.yaml")
    preflight.to_csv(temporary / "preflight.csv", index=False)
    (temporary / "setup_report.json").write_text(
        json.dumps(dict(setup_report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_fasta(temporary / "reference_sequences.fasta", reference_sequences)
    for source in structure_paths:
        target = temporary / "structures" / source.name
        if target.exists() and _sha256(target) != _sha256(source):
            raise ValueError(f"Structure filename collision with different content: {source.name}")
        shutil.copy2(source, target)
        manifest["structures"].append(f"structures/{source.name}")
    for path in sorted(temporary.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            relative = str(path.relative_to(temporary))
            manifest["files"][relative] = {
                "sha256": _sha256(path),
                "size": path.stat().st_size,
            }
    (temporary / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if destination.exists():
        destination.unlink()
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(temporary.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=str(path.relative_to(temporary)))
    shutil.rmtree(temporary)
    return destination


def extract_and_verify_setup_bundle(bundle: str | Path, destination: str | Path) -> dict[str, Any]:
    bundle = Path(bundle)
    destination = Path(destination)
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    with zipfile.ZipFile(bundle) as archive:
        members = archive.namelist()
        unsafe = [name for name in members if Path(name).is_absolute() or ".." in Path(name).parts]
        if unsafe:
            raise ValueError(f"Unsafe bundle paths: {unsafe}")
        archive.extractall(destination)
    manifest_path = destination / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("Setup bundle lacks manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("setup_wizard_version") != SETUP_WIZARD_VERSION:
        raise ValueError(
            "Setup bundle version differs: " + repr(manifest.get("setup_wizard_version"))
        )
    for relative, record in manifest.get("files", {}).items():
        path = destination / relative
        if not path.is_file():
            raise ValueError(f"Setup bundle file is missing: {relative}")
        if _sha256(path) != record.get("sha256"):
            raise ValueError(f"Setup bundle hash mismatch: {relative}")
        if path.stat().st_size != int(record.get("size", -1)):
            raise ValueError(f"Setup bundle size mismatch: {relative}")
    return manifest


def _systems_mapping(config: Mapping[str, Any]) -> dict[str, Any]:
    systems = config.get("systems")
    if isinstance(systems, dict):
        return dict(systems)
    if isinstance(systems, list):
        result = {}
        for item in systems:
            if not isinstance(item, dict) or "name" not in item:
                raise ValueError("List-form config system lacks a name.")
            name = str(item["name"])
            value = dict(item)
            value.pop("name", None)
            result[name] = value
        return result
    raise ValueError("Config lacks a systems mapping or list.")


def merge_setup_bundles(
    bundles: Sequence[str | Path],
    destination: str | Path,
) -> tuple[Path, Path, dict[str, Any]]:
    """Merge one or more validated setup bundles for multi-system batch runs."""
    destination = Path(destination)
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    structures_dir = destination / "structures"
    structures_dir.mkdir()
    combined_systems: dict[str, Any] = {}
    combined_sequences: dict[str, str] = {}
    source_records: list[dict[str, Any]] = []
    for index, bundle in enumerate(bundles, start=1):
        extracted = destination / f"bundle_{index}"
        manifest = extract_and_verify_setup_bundle(bundle, extracted)
        config = yaml.safe_load((extracted / manifest["config"]).read_text(encoding="utf-8"))
        systems = _systems_mapping(config)
        bundle_sequences = read_fasta(extracted / "reference_sequences.fasta")
        for gene, sequence in bundle_sequences.items():
            if gene in combined_sequences and combined_sequences[gene] != sequence:
                raise ValueError(
                    f"Prepared bundles contain conflicting reference sequences for {gene!r}."
                )
            combined_sequences[gene] = sequence
        overlap = sorted(set(combined_systems) & set(systems))
        if overlap:
            raise ValueError(f"Duplicate system names across bundles: {overlap}")
        combined_systems.update(systems)
        for relative in manifest.get("structures", []):
            source = extracted / relative
            target = structures_dir / source.name
            if target.exists() and _sha256(target) != _sha256(source):
                raise ValueError(
                    f"Structure filename collision with different content: {source.name}"
                )
            shutil.copy2(source, target)
        source_records.append(
            {
                "bundle": Path(bundle).name,
                "systems": sorted(systems),
                "manifest_sha256": _sha256(extracted / "manifest.json"),
            }
        )
    config_path = destination / "config.yaml"
    config_path.write_text(
        yaml.safe_dump({"systems": combined_systems}, sort_keys=False),
        encoding="utf-8",
    )
    reference_fasta = write_fasta(
        destination / "reference_sequences.fasta",
        combined_sequences,
    )
    merge_manifest = {
        "setup_wizard_version": SETUP_WIZARD_VERSION,
        "systems": sorted(combined_systems),
        "sources": source_records,
        "reference_sequences": reference_fasta.name,
        "genes": sorted(combined_sequences),
    }
    (destination / "merge_manifest.json").write_text(
        json.dumps(merge_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return config_path, structures_dir, merge_manifest


def _system_gene_mapping(system: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    genes = system.get("genes", {})
    if isinstance(genes, dict):
        return {normalize_gene(gene): dict(value) for gene, value in genes.items()}
    if isinstance(genes, list):
        result: dict[str, dict[str, Any]] = {}
        for item in genes:
            if not isinstance(item, dict) or not item.get("gene"):
                raise ValueError("List-form gene configuration lacks a gene identifier.")
            value = dict(item)
            gene = normalize_gene(value.pop("gene"))
            result[gene] = value
        return result
    raise ValueError("System gene configuration is neither a mapping nor a list.")


def _system_complex_filename(system: Mapping[str, Any]) -> str:
    for key in ("pdb_file", "cif_file", "complex_file", "structure_file"):
        value = system.get(key)
        if value:
            return Path(str(value)).name
    raise ValueError("System configuration does not identify a complex coordinate file.")


def _gene_monomer_filename(gene_config: Mapping[str, Any]) -> str:
    for key in ("monomer_pdb", "monomer_cif", "monomer_file"):
        value = gene_config.get(key)
        if value:
            return Path(str(value)).name
    raise ValueError("Gene configuration does not identify a monomer coordinate file.")


def revalidate_prepared_systems(
    config_path: str | Path,
    structures_dir: str | Path,
    reference_sequences: Mapping[str, str],
    variants: pd.DataFrame,
    *,
    output_config_path: str | Path,
    input_offset_overrides: Mapping[str, int] | None = None,
) -> tuple[Path, pd.DataFrame, dict[str, Any]]:
    """Recheck prepared systems against the current variant batch.

    A prepared bundle stores reviewed structures and chain identities, but its
    combined COMAVI offsets reflect the numbering convention used when it was
    created.  This function keeps the reviewed chains and structures while
    re-inferring the submitted-numbering convention and rebuilding the offsets
    for the current batch.  Every current variant is rechecked against both the
    monomer and complex before FoldX.
    """
    config_path = Path(config_path)
    structures_dir = Path(structures_dir)
    output_config_path = Path(output_config_path)
    frame = _variant_frame(variants)
    references = {
        normalize_gene(gene): normalize_sequence(sequence)
        for gene, sequence in reference_sequences.items()
    }
    input_offset_overrides = {
        normalize_gene(gene): int(value)
        for gene, value in (input_offset_overrides or {}).items()
    }
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    systems = _systems_mapping(config)
    updated_systems: dict[str, Any] = {}
    preflight_frames: list[pd.DataFrame] = []
    reports: dict[str, Any] = {}

    for system_name, raw_system in systems.items():
        system = dict(raw_system)
        gene_config = _system_gene_mapping(system)
        system_genes = list(gene_config)
        missing_references = [gene for gene in system_genes if gene not in references]
        if missing_references:
            raise ValueError(
                f"Prepared system {system_name!r} lacks reference sequences for: "
                + ", ".join(missing_references)
            )
        system_variants = frame.loc[frame["gene"].isin(system_genes)].copy()
        if system_variants.empty:
            updated_systems[str(system_name)] = system
            reports[str(system_name)] = {
                "status": "not_used",
                "message": "No submitted variant belongs to this prepared system.",
            }
            continue

        complex_path = structures_dir / _system_complex_filename(system)
        chain_overrides = {
            gene: str(settings.get("chain"))
            for gene, settings in gene_config.items()
            if settings.get("chain") not in (None, "")
        }
        complex_assessment = assess_structure(
            complex_path,
            {gene: references[gene] for gene in system_genes},
            system_variants,
            source_label=f"prepared:{system_name}",
            input_offset_overrides=input_offset_overrides,
            chain_overrides=chain_overrides,
        )

        monomer_assessments: dict[str, MonomerAssessment] = {}
        for gene, settings in gene_config.items():
            monomer_path = structures_dir / _gene_monomer_filename(settings)
            monomer_assessments[gene] = assess_monomer(
                gene,
                monomer_path,
                references[gene],
                system_variants,
                input_offset_override=input_offset_overrides.get(gene),
            )

        current_variant_genes = set(system_variants["gene"])
        updated_gene_config: dict[str, dict[str, Any]] = {}
        for gene, settings in gene_config.items():
            updated = dict(settings)
            if gene in current_variant_genes:
                monomer_offset = monomer_assessments[gene].pipeline_offset
                multimer_offset = complex_assessment.pipeline_offsets.get(gene)
                if monomer_offset is None or multimer_offset is None:
                    raise ValueError(
                        f"Prepared system {system_name!r} could not derive current offsets for {gene}."
                    )
                updated["monomer_offset"] = int(monomer_offset)
                updated["multimer_offset"] = int(multimer_offset)
            updated_gene_config[gene] = updated
        system["genes"] = updated_gene_config
        updated_systems[str(system_name)] = system

        preflight = build_preflight_table(
            system_variants,
            monomer_assessments,
            complex_assessment,
        )
        preflight.insert(0, "system", str(system_name))
        preflight_frames.append(preflight)
        reports[str(system_name)] = {
            "status": (
                "red"
                if (preflight["overall"] == "red").any()
                else "yellow"
                if (preflight["overall"] == "yellow").any()
                else "green"
            ),
            "complex_assessment": assessment_to_jsonable(complex_assessment),
            "monomer_assessments": {
                gene: assessment_to_jsonable(value)
                for gene, value in monomer_assessments.items()
            },
        }

    output_config_path.parent.mkdir(parents=True, exist_ok=True)
    output_config_path.write_text(
        yaml.safe_dump({"systems": updated_systems}, sort_keys=False),
        encoding="utf-8",
    )
    if preflight_frames:
        combined_preflight = pd.concat(preflight_frames, ignore_index=True)
    else:
        combined_preflight = pd.DataFrame(
            columns=[
                "system", "gene", "variant", "input", "monomer_mapping",
                "complex_mapping", "overall", "action",
            ]
        )
    return output_config_path, combined_preflight, reports


def genes_from_config(config_path: str | Path) -> list[str]:
    config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    systems = _systems_mapping(config)
    genes: list[str] = []
    for system in systems.values():
        gene_block = system.get("genes", {}) if isinstance(system, dict) else {}
        if isinstance(gene_block, dict):
            candidates = gene_block.keys()
        elif isinstance(gene_block, list):
            candidates = [item.get("gene") for item in gene_block if isinstance(item, dict)]
        else:
            candidates = []
        for gene in candidates:
            key = normalize_gene(gene)
            if key and key not in genes:
                genes.append(key)
    return genes


def assessment_to_jsonable(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        return {key: assessment_to_jsonable(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): assessment_to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [assessment_to_jsonable(item) for item in value]
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    return value
