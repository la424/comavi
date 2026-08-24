#!/usr/bin/env python3
"""Verify the CHD structure manifest against the config and optional files."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import yaml


REPO = Path(__file__).resolve().parents[1]

DEFAULT_CONFIG = (
    REPO
    / "configs"
    / "chd_systems.yaml"
)

DEFAULT_MANIFEST = (
    REPO
    / "reference_outputs"
    / "chd_isds_v1"
    / "chd_structure_manifest.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="CHD systems YAML configuration.",
    )

    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Path-free CHD structure manifest.",
    )

    parser.add_argument(
        "--structure-dir",
        type=Path,
        help=(
            "Optional directory containing the actual "
            "structure files. When supplied, sizes and "
            "SHA-256 hashes are verified."
        ),
    )

    return parser.parse_args()


def require(
    condition: bool,
    message: str,
) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def sha256(file_path: Path) -> str:
    digest = hashlib.sha256()

    with file_path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def derive_requirements(
    config_value: dict[str, Any],
) -> dict[str, dict[str, set[str]]]:
    systems = config_value.get(
        "systems",
        {},
    )

    require(
        isinstance(systems, dict)
        and bool(systems),
        "CHD configuration lacks systems.",
    )

    requirements: dict[
        str,
        dict[str, set[str]],
    ] = defaultdict(
        lambda: {
            "kinds": set(),
            "roles": set(),
            "systems": set(),
        }
    )

    def add(
        filename: Any,
        *,
        kind: str,
        role: str,
        system: str,
    ) -> None:
        if filename is None:
            return

        text = str(filename).strip()

        if not text:
            return

        require(
            Path(text).name == text,
            f"configured structure is not a basename: {text}",
        )

        record = requirements[text]
        record["kinds"].add(kind)
        record["roles"].add(role)
        record["systems"].add(system)

    for system, specification in systems.items():
        require(
            isinstance(specification, dict),
            f"invalid system specification: {system}",
        )

        add(
            specification.get("pdb_file"),
            kind="PDB",
            role="multimer",
            system=str(system),
        )

        add(
            specification.get("cif_file"),
            kind="CIF",
            role="multimer-confidence",
            system=str(system),
        )

        genes = specification.get(
            "genes",
            {},
        )

        require(
            isinstance(genes, dict)
            and bool(genes),
            f"{system} lacks gene definitions.",
        )

        for gene, gene_specification in genes.items():
            require(
                isinstance(
                    gene_specification,
                    dict,
                ),
                f"invalid gene specification: {system}/{gene}",
            )

            add(
                gene_specification.get(
                    "monomer_pdb"
                ),
                kind="PDB",
                role=f"monomer:{gene}",
                system=str(system),
            )

            add(
                gene_specification.get(
                    "monomer_cif"
                ),
                kind="CIF",
                role=(
                    f"monomer-confidence:{gene}"
                ),
                system=str(system),
            )

    return dict(requirements)


def main() -> None:
    args = parse_args()

    config_file = args.config.expanduser().resolve()
    manifest_file = args.manifest.expanduser().resolve()

    require(
        config_file.is_file(),
        f"configuration is absent: {config_file}",
    )

    require(
        manifest_file.is_file(),
        f"manifest is absent: {manifest_file}",
    )

    config_text = config_file.read_text(
        encoding="utf-8"
    )

    config_value = yaml.safe_load(
        config_text
    )

    require(
        isinstance(config_value, dict),
        "configuration did not parse as a mapping.",
    )

    requirements = derive_requirements(
        config_value
    )

    manifest = json.loads(
        manifest_file.read_text(
            encoding="utf-8"
        )
    )

    require(
        isinstance(manifest, dict),
        "structure manifest is not an object.",
    )

    require(
        manifest.get("schema_version")
        == "CHD-structure-manifest-v1",
        "unexpected structure-manifest version.",
    )

    expected_config_hash = hashlib.sha256(
        config_text.encode("utf-8")
    ).hexdigest()

    require(
        manifest.get("config_sha256")
        == expected_config_hash,
        "manifest config SHA-256 differs from the current config.",
    )

    authority_commit = str(
        manifest.get(
            "authority_commit",
            "",
        )
    )

    require(
        re.fullmatch(
            r"[0-9a-f]{40}",
            authority_commit,
        )
        is not None,
        "invalid authority commit.",
    )

    entries = manifest.get("files")

    require(
        isinstance(entries, list),
        "manifest files field is not a list.",
    )

    require(
        int(
            manifest.get(
                "file_count",
                -1,
            )
        )
        == len(entries),
        "manifest file_count differs from the entry count.",
    )

    manifest_by_name: dict[
        str,
        dict[str, Any],
    ] = {}

    for entry in entries:
        require(
            isinstance(entry, dict),
            "a manifest file entry is not an object.",
        )

        require(
            set(entry)
            == {
                "filename",
                "size_bytes",
                "sha256",
            },
            "a manifest entry has unexpected fields.",
        )

        filename = str(
            entry["filename"]
        ).strip()

        require(
            filename
            and Path(filename).name
            == filename,
            f"invalid path-free filename: {filename!r}",
        )

        key = filename.casefold()

        require(
            key not in manifest_by_name,
            f"duplicate manifest filename: {filename}",
        )

        size_bytes = entry["size_bytes"]
        digest = str(
            entry["sha256"]
        ).lower()

        require(
            isinstance(size_bytes, int)
            and size_bytes > 0,
            f"invalid file size for {filename}",
        )

        require(
            re.fullmatch(
                r"[0-9a-f]{64}",
                digest,
            )
            is not None,
            f"invalid SHA-256 for {filename}",
        )

        manifest_by_name[key] = {
            "filename": filename,
            "size_bytes": size_bytes,
            "sha256": digest,
        }

    required_by_name = {
        filename.casefold(): filename
        for filename in requirements
    }

    require(
        set(manifest_by_name)
        == set(required_by_name),
        (
            "manifest filenames differ from the current "
            "configuration requirements.\n"
            f"missing={sorted(set(required_by_name) - set(manifest_by_name))}\n"
            f"extra={sorted(set(manifest_by_name) - set(required_by_name))}"
        ),
    )

    pdb_count = sum(
        "PDB" in record["kinds"]
        for record in requirements.values()
    )

    cif_count = sum(
        "CIF" in record["kinds"]
        for record in requirements.values()
    )

    if args.structure_dir is not None:
        structure_dir = (
            args.structure_dir
            .expanduser()
            .resolve()
        )

        require(
            structure_dir.is_dir(),
            f"structure directory is absent: {structure_dir}",
        )

        for key in sorted(manifest_by_name):
            entry = manifest_by_name[key]
            structure_file = (
                structure_dir
                / entry["filename"]
            )

            require(
                structure_file.is_file(),
                f"required structure is absent: {entry['filename']}",
            )

            require(
                structure_file.stat().st_size
                == entry["size_bytes"],
                f"size differs for {entry['filename']}",
            )

            require(
                sha256(structure_file)
                == entry["sha256"],
                f"SHA-256 differs for {entry['filename']}",
            )

    print(f"Configured systems: {len(config_value['systems'])}")
    print(f"Required structure files: {len(requirements)}")
    print(f"Required PDB files: {pdb_count}")
    print(f"Required CIF files: {cif_count}")
    print(f"Authority commit: {authority_commit}")

    if args.structure_dir is not None:
        print(
            "External structure package: "
            "all sizes and hashes match"
        )
    else:
        print(
            "Static config-manifest contract: match"
        )

    print("CHD STRUCTURE-CONTRACT VERIFICATION: PASS")


if __name__ == "__main__":
    main()
