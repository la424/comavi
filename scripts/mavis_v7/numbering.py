"""
Residue-identity check for offset calibration.

Why this module exists
----------------------
Every variant is scored at two independent structure positions:

    position_mono  = CSV position - monomers[gene].position_offset
    position_multi = CSV position - multimer.position_offsets[gene]

Those two offsets are calibrated against *different* coordinate sources — the
monomer axis against a predicted model (AlphaFold numbering, which retains the
initiator Met), the complex axis against a crystal or NMR entry (author
numbering, which often does not, and which may start partway into the protein
for a domain construct). Getting one wrong does not raise: FoldX happily mutates
whatever residue sits at the position it is handed.

That failure is silent whenever the wrong position holds the same amino acid as
the intended one, which is not rare. The worked example in this repo is exactly
that case: HBB E6V with `monomer_offset: -1` resolves to AlphaFold residue 7 =
Glu, and crystal chain B residue 7 is *also* Glu, so using crystal chains for
the monomer axis scores Glu7 instead of Glu6 and returns a plausible number.

So the reference amino acid in the variant CSV is a checksum on the offsets, and
this module spends it. Run before FoldX, not after.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .config import SystemConfig
from .structure_loading import load_pdb, get_residue_aa

# How far either side of the configured offset to look when suggesting a fix.
OFFSET_SEARCH = 25


class NumberingIssue:
    """One (variant, axis) mismatch, with the offsets that would have worked."""

    def __init__(self, system, gene, variant, axis, path, chain,
                 position_csv, position_struct, expected_aa, found_aa,
                 configured_offset, candidate_offsets, silent_trap):
        self.system = system
        self.gene = gene
        self.variant = variant
        self.axis = axis                      # "monomer" | "complex"
        self.path = path
        self.chain = chain
        self.position_csv = position_csv
        self.position_struct = position_struct
        self.expected_aa = expected_aa
        self.found_aa = found_aa              # None => position absent from chain
        self.configured_offset = configured_offset
        self.candidate_offsets = candidate_offsets
        self.silent_trap = silent_trap        # found_aa == expected_aa at wrong position

    def __str__(self):
        key = "monomer_offset" if self.axis == "monomer" else "multimer_offset"
        if self.found_aa is None:
            what = (f"residue {self.position_struct} is absent from chain "
                    f"{self.chain} (construct may not span this position)")
        else:
            what = (f"residue {self.position_struct} in chain {self.chain} is "
                    f"{self.found_aa}, but the variant says {self.expected_aa}")
        lines = [f"[{self.system}] {self.gene} {self.variant} - {self.axis} axis: {what}",
                 f"    structure        {self.path}",
                 f"    CSV position     {self.position_csv}",
                 f"    {key:16s} {self.configured_offset}  "
                 f"-> structure position {self.position_struct}"]
        if self.candidate_offsets:
            shown = ", ".join(str(o) for o in self.candidate_offsets[:6])
            lines.append(f"    offsets that would put {self.expected_aa} here: {shown}")
        else:
            lines.append(f"    no offset within +/-{OFFSET_SEARCH} places "
                         f"{self.expected_aa} at this position - wrong chain or "
                         f"wrong structure?")
        return "\n".join(lines)


def _candidate_offsets(aa_map: Dict[int, str], position_csv: int,
                       expected_aa: str) -> List[int]:
    """Offsets o (position_csv - o present in aa_map as expected_aa), nearest first."""
    hits = [o for o in range(-OFFSET_SEARCH, OFFSET_SEARCH + 1)
            if aa_map.get(position_csv - o) == expected_aa]
    return sorted(hits, key=abs)


def _resolve_multimer_path(cfg: SystemConfig, structure_dir: Path,
                           preprocessed_dir: Path) -> Optional[Path]:
    """Preprocessed complex if it exists, else the raw source.

    Preprocessing strips HETATM records but does not renumber, so the raw source
    is a valid substitute for a numbering check before preprocessing has run.
    """
    multi = cfg.multimer
    if multi.structure_type == "AF":
        p = Path(structure_dir) / multi.pdb_file
        return p if p.exists() else None
    processed = Path(preprocessed_dir) / multi.pdb_file
    if processed.exists():
        return processed
    raw = Path(structure_dir) / (multi.source_file or multi.pdb_file)
    return raw if raw.exists() else None


def _check_one(aa_map, system, gene, variant, axis, path, chain,
               position_csv, position_struct, expected_aa, configured_offset):
    found = aa_map.get(position_struct)
    if found == expected_aa:
        return None
    cands = _candidate_offsets(aa_map, position_csv, expected_aa)
    return NumberingIssue(
        system=system, gene=gene, variant=variant, axis=axis, path=path,
        chain=chain, position_csv=position_csv, position_struct=position_struct,
        expected_aa=expected_aa, found_aa=found,
        configured_offset=configured_offset, candidate_offsets=cands,
        silent_trap=bool(cands),
    )


def check_variant_numbering(
    configs: Dict[str, SystemConfig],
    variants,
    structure_dir,
    preprocessed_dir,
) -> Tuple[bool, List[NumberingIssue], Dict[str, int]]:
    """Verify the reference AA sits where each axis's offset says it does.

    `variants` is the expanded frame: it must carry gene, system, ref_aa,
    position, position_mono, position_multi (see variant_loading.py).

    Returns (all_ok, issues, counts) where counts summarises what was checked.
    Structures that cannot be loaded are counted as skipped, not as failures -
    file existence is validate_config's job, not this function's.
    """
    structure_dir = Path(structure_dir)
    preprocessed_dir = Path(preprocessed_dir)
    issues: List[NumberingIssue] = []
    counts = {"checked": 0, "skipped": 0, "monomer": 0, "complex": 0,
              "unresolvable": 0}
    cache: Dict[Tuple[str, str], Dict[int, str]] = {}

    def aa_map_for(path: Path, chain: str) -> Optional[Dict[int, str]]:
        key = (str(path), chain)
        if key not in cache:
            struct = load_pdb(path)
            cache[key] = get_residue_aa(struct, chain) if struct is not None else {}
        return cache[key] or None

    for _, row in variants.iterrows():
        system, gene = row["system"], row["gene"]
        cfg = configs.get(system)
        if cfg is None:
            continue
        variant = f"{row['ref_aa']}{row['position']}{row['alt_aa']}"
        expected = str(row["ref_aa"]).strip().upper()

        # ---- monomer axis ----
        mspec = cfg.monomers.get(gene)
        if mspec is None:
            counts["unresolvable"] += 1
        else:
            path = structure_dir / mspec.pdb_file
            aa_map = aa_map_for(path, mspec.chain_id)
            if aa_map is None:
                counts["skipped"] += 1
            else:
                counts["checked"] += 1
                counts["monomer"] += 1
                iss = _check_one(
                    aa_map, system, gene, variant, "monomer", path,
                    mspec.chain_id, int(row["position"]),
                    int(row["position_mono"]), expected, mspec.position_offset)
                if iss:
                    issues.append(iss)

        # ---- complex axis ----
        chain = cfg.multimer.chain_map.get(gene)
        path = _resolve_multimer_path(cfg, structure_dir, preprocessed_dir)
        if chain and path:
            aa_map = aa_map_for(path, chain)
            if aa_map is None:
                counts["skipped"] += 1
            else:
                counts["checked"] += 1
                counts["complex"] += 1
                iss = _check_one(
                    aa_map, system, gene, variant, "complex", path, chain,
                    int(row["position"]), int(row["position_multi"]), expected,
                    cfg.multimer.position_offsets.get(gene, 0))
                if iss:
                    issues.append(iss)
        else:
            counts["skipped"] += 1

    # "No mismatches" only means something if positions were actually verified.
    ok = len(issues) == 0 and counts["checked"] > 0
    return ok, issues, counts


def format_report(ok: bool, issues: List[NumberingIssue],
                  counts: Dict[str, int]) -> str:
    """Human-readable block for run.py / scripts to print."""
    if ok and counts["checked"] == 0:
        # A check that verified nothing must not report success -- that is a false
        # green light, and it is what a gene-name or path mismatch looks like.
        return ("Residue-identity check INCONCLUSIVE: 0 axis positions were "
                f"verified ({counts['skipped']} skipped, "
                f"{counts.get('unresolvable', 0)} gene(s) not resolvable against the "
                "config). Usually a gene-name mismatch between the variant frame and "
                "the config's gene keys, or structures that did not load.")
    if ok:
        return (f"Residue-identity check passed: {counts['checked']} axis position(s) "
                f"({counts['monomer']} monomer, {counts['complex']} complex)"
                + (f", {counts['skipped']} skipped (structure not loadable)"
                   if counts["skipped"] else "")
                + (f", {counts['unresolvable']} gene(s) unresolvable"
                   if counts.get("unresolvable") else ""))
    lines = [f"Residue-identity check FAILED: {len(issues)} mismatch(es) "
             f"across {counts['checked']} axis position(s) checked", ""]
    for iss in issues:
        lines.append(str(iss))
        lines.append("")
    if any(i.silent_trap for i in issues):
        lines.append(
            "At least one mismatch would be SILENT: another offset places the "
            "expected residue nearby, so FoldX would have returned a "
            "plausible-looking value for the wrong residue. Fix the offsets "
            "before scoring.")
    return "\n".join(lines)


# ============================================================================
# Structure provenance
# ============================================================================
# The residue-identity check above cannot catch every offset error, and the
# repo's own worked example is the proof: HBB E6V with monomer_offset -1 resolves
# to residue 7, and residue 7 is Glu in BOTH the AlphaFold model (correct) and
# crystal chain B (wrong). Substituting crystal chains for the monomer models
# passes the identity check and still scores the wrong residue.
#
# What separates the two sources is the B-factor column. In a predicted model it
# carries pLDDT, which is a per-residue quantity, so every atom in a residue holds
# an identical value. In an experimental structure it carries a real per-atom
# temperature factor, which essentially never repeats exactly across a residue.
# Measured on this repo's structures: AlphaFold models 100% uniform, crystal
# entries 0%.
#
# This check has a second payoff independent of numbering. pLDDT gating reads the
# same column, so an experimental file standing in for a predicted monomer would
# have its crystallographic B-factors interpreted as confidence scores.

UNIFORM_B_PREDICTED = 0.90   # >= this fraction of residues uniform => predicted
UNIFORM_B_EXPERIMENTAL = 0.10  # <= this fraction => experimental
_B_TOL = 1e-6


def b_factor_uniformity(path) -> Tuple[Optional[float], int]:
    """Fraction of residues whose atoms share one B-factor, and residue count.

    Returns (None, 0) if the file has no parsable ATOM records. Single-atom
    residues are excluded: they are trivially uniform and would bias the ratio.
    """
    path = Path(path)
    if not path.exists():
        return None, 0
    per_res: Dict[Tuple[str, int, str], List[float]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            b = float(line[60:66])
        except (ValueError, IndexError):
            continue
        per_res.setdefault((line[21], int(line[22:26]), line[26]), []).append(b)
    multi = [v for v in per_res.values() if len(v) > 1]
    if not multi:
        return None, len(per_res)
    uniform = sum(1 for v in multi if (max(v) - min(v)) < _B_TOL)
    return uniform / len(multi), len(per_res)


def classify_structure_source(path) -> Optional[str]:
    """'predicted' | 'experimental' | 'uninformative' | 'ambiguous'.

    None if the file has no parsable ATOM records.

    'uninformative' matters and was found the hard way: NMR entries frequently
    carry an all-zero B-factor column (1JM7, this benchmark's BRCA1-BARD1
    structure, is one). A single constant value across the whole file is
    per-residue-uniform by construction, so a naive uniformity test calls it
    'predicted'. It is not evidence of anything -- the column is simply empty --
    and treating it as predicted would fail a correctly configured NMR system.
    A real predicted model carries pLDDT, which varies from residue to residue.
    """
    frac, _ = b_factor_uniformity(path)
    if frac is None:
        return None
    p = Path(path)
    distinct = set()
    for line in p.read_text().splitlines():
        if line.startswith("ATOM"):
            b = line[60:66].strip()
            if b:
                distinct.add(b)
                if len(distinct) > 1:
                    break
    if len(distinct) <= 1:
        return "uninformative"
    if frac >= UNIFORM_B_PREDICTED:
        return "predicted"
    if frac <= UNIFORM_B_EXPERIMENTAL:
        return "experimental"
    return "ambiguous"


def check_structure_provenance(configs: Dict[str, SystemConfig],
                               structure_dir) -> Tuple[bool, List[str]]:
    """Confirm each monomer file's source matches what its spec declares.

    Returns (all_ok, messages). A monomer spec with structure_type 'AF' (the
    default) or plddt_gate on is expected to point at a predicted model; a file
    that looks experimental is an error, because it means both the numbering
    offset and the pLDDT gate are being applied to the wrong kind of coordinates.
    """
    structure_dir = Path(structure_dir)
    msgs: List[str] = []
    ok = True
    seen: Dict[str, Optional[str]] = {}

    for sys_name, cfg in configs.items():
        for gene, mspec in cfg.monomers.items():
            path = structure_dir / mspec.pdb_file
            if str(path) not in seen:
                seen[str(path)] = classify_structure_source(path)
            kind = seen[str(path)]
            if kind is None:
                continue
            wants_predicted = (mspec.structure_type == "AF") or mspec.plddt_gate
            if wants_predicted and kind == "experimental":
                ok = False
                msgs.append(
                    f"[{sys_name}] monomer {gene}: {path.name} looks EXPERIMENTAL "
                    f"(per-atom B-factors) but the spec declares "
                    f"structure_type={mspec.structure_type!r}, plddt_gate="
                    f"{mspec.plddt_gate}, i.e. a predicted model.\n"
                    f"    Two things break: monomer_offset="
                    f"{mspec.position_offset} is calibrated to predicted-model "
                    f"numbering, and pLDDT gating will read crystallographic "
                    f"B-factors as confidence scores.\n"
                    f"    If you meant to score a monomer taken from a crystal, set "
                    f"structure_type: xray and plddt_gate: false on this gene, and "
                    f"recalibrate monomer_offset against that entry's numbering.")
            elif not wants_predicted and kind == "predicted":
                msgs.append(
                    f"[{sys_name}] monomer {gene}: {path.name} looks PREDICTED "
                    f"(uniform per-residue B-factors) but the spec declares "
                    f"structure_type={mspec.structure_type!r}. Check "
                    f"monomer_offset is calibrated to this file.")

        # ---- multimer: same logic, and the same pLDDT-gate consequence ----
        multi = getattr(cfg, "multimer", None)
        if multi is None:
            continue
        mpath = structure_dir / multi.pdb_file
        if str(mpath) not in seen:
            seen[str(mpath)] = classify_structure_source(mpath)
        kind = seen[str(mpath)]
        if kind is None:
            continue
        wants_predicted = (multi.structure_type == "AF") or multi.plddt_gate
        if wants_predicted and kind == "experimental":
            ok = False
            msgs.append(
                f"[{sys_name}] multimer: {mpath.name} looks EXPERIMENTAL but the "
                f"spec declares structure_type={multi.structure_type!r}, "
                f"plddt_gate={multi.plddt_gate}. Interface pLDDT gating would read "
                f"crystallographic B-factors as confidence scores. Set "
                f"structure_type: xray and plddt_gate: false for this system.")
        elif not wants_predicted and kind == "predicted":
            # Warn, do not fail. Scoring a predicted interface with gating off is a
            # deliberate choice in some configs, and the direction that silently
            # corrupts results is the other one (experimental file where a
            # predicted model is expected).
            msgs.append(
                f"[{sys_name}] multimer: {mpath.name} looks PREDICTED (uniform "
                f"per-residue B-factors) but the spec declares "
                f"structure_type={multi.structure_type!r} with plddt_gate="
                f"{multi.plddt_gate}. A predicted interface is being scored with "
                f"confidence gating disabled, and multimer_offset may be "
                f"calibrated to the wrong numbering.")
    return ok, msgs
