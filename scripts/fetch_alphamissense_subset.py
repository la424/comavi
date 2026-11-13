#!/usr/bin/env python3
"""Extract the benchmark proteins from the AlphaMissense data release.

WHY THIS EXISTS
---------------
The 44 AlphaMissense scores in inputs/AM_variants_comavi_mechanism_test.xlsx
and the three CFH scores in the canonical table were entered without a
recorded source, and the two VWF variants (p.Arg1334Gln, p.Thr1381Ala) had
none at all. AlphaFold DB serves per-protein AlphaMissense tables only for
proteins it hosts, and holds only isoform 2 of VWF (residues 1-351), so the
scores have to come from the release file itself.

This script streams the published release once, verifies the MD5 published
with the record, and keeps the rows for the benchmark proteins. The subset is
committed so that every AlphaMissense value in the benchmark can be traced to a
checksummed public file without re-downloading 1.2 GB.

Source: Cheng J et al. (2023) Science 381:eadg7492. Predictions for
AlphaMissense, Zenodo record 8208688 (doi:10.5281/zenodo.8208688), file
AlphaMissense_aa_substitutions.tsv.gz. Coordinates are UniProt canonical
positions counted from the initiator methionine.

Usage:
    python scripts/fetch_alphamissense_subset.py
"""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib
import sys
import time
import urllib.request
import zlib

REPO = pathlib.Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "inputs" / "alphamissense"
URL = ("https://zenodo.org/records/8208688/files/"
       "AlphaMissense_aa_substitutions.tsv.gz")
EXPECTED_MD5 = "b9ccb339e0de6cb0a8d1973ad2026576"
EXPECTED_BYTES = 1207278510

# UniProt canonical accessions (reviewed, human) for every benchmark gene that
# carries a clinical label. Resolved with rest.uniprot.org gene_exact queries.
ACCESSIONS = {
    "BRCA1": "P38398", "CALM1": "P0DP23", "CFH": "P08603", "HBB": "P68871",
    "KRAS": "P01116", "MLH1": "P40692", "MSH2": "P43246", "PIK3CA": "P42336",
    "PIK3R1": "P27986", "PMS2": "P54278", "SMAD4": "Q13485", "TNNI3": "P19429",
    "VHL": "P40337", "VWF": "P04275",
}


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    want = {a.encode() for a in ACCESSIONS.values()}
    out_path = OUT_DIR / "AlphaMissense_benchmark_subset.tsv"
    md5 = hashlib.md5()
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    tail = b""
    n_bytes = kept = 0
    header = None
    t0 = time.time()
    next_report = 100 << 20
    with urllib.request.urlopen(URL, timeout=180) as resp, out_path.open("w") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            md5.update(chunk)
            n_bytes += len(chunk)
            data = dec.decompress(chunk)
            # The release may be multi-member gzip: restart on each member.
            while dec.eof and dec.unused_data:
                rest = dec.unused_data
                dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
                data += dec.decompress(rest)
            lines = (tail + data).split(b"\n")
            tail = lines.pop()
            for ln in lines:
                if ln.startswith(b"#") or ln.startswith(b"uniprot_id"):
                    if b"uniprot_id" in ln:
                        header = ln.lstrip(b"#").decode()
                        out.write(header + "\n")
                    continue
                if ln.split(b"\t", 1)[0] in want:
                    out.write(ln.decode() + "\n")
                    kept += 1
            if n_bytes >= next_report:
                rate = n_bytes / max(time.time() - t0, 1e-6) / (1 << 20)
                print("  %5d MB of %d  (%.1f MB/s)  rows kept %d"
                      % (n_bytes >> 20, EXPECTED_BYTES >> 20, rate, kept), flush=True)
                next_report += 100 << 20
        if tail and tail.split(b"\t", 1)[0] in want:
            out.write(tail.decode() + "\n")
            kept += 1

    digest = md5.hexdigest()
    ok = digest == EXPECTED_MD5 and n_bytes == EXPECTED_BYTES
    prov = {
        "source": "Predictions for AlphaMissense, Zenodo record 8208688",
        "doi": "10.5281/zenodo.8208688",
        "citation": "Cheng J et al. (2023) Science 381:eadg7492",
        "file": "AlphaMissense_aa_substitutions.tsv.gz",
        "url": URL,
        "bytes": n_bytes,
        "md5": digest,
        "md5_expected": EXPECTED_MD5,
        "md5_verified": ok,
        "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "accessions": ACCESSIONS,
        "header": header,
        "rows_kept": kept,
        "coordinates": "UniProt canonical positions from the initiator methionine",
    }
    (OUT_DIR / "AlphaMissense_benchmark_subset.provenance.json").write_text(
        json.dumps(prov, indent=1) + "\n")
    print("bytes %d | md5 %s | verified %s | rows kept %d | %.0f s"
          % (n_bytes, digest, ok, kept, time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
