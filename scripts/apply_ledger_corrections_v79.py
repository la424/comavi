#!/usr/bin/env python3
"""v7.9 -- Reverse three v7.6 token changes whose basis contradicts the token.

WHY THIS EXISTS
---------------
An audit of every committed expectation against its own written evidence basis
found three rows whose `expected_token` asserts the opposite of the basis
recorded beside it:

  msh2_msh6 C697F  binding  destab  <- basis: "ATPase catalytic defect,
                                       protein stable, interaction intact"
  kras_craf G12D   binding  destab  <- basis: "WT-like RAF1 binding; impaired
                                       GTP hydrolysis. G12 is 13.6 A from RAF1"
  kras_craf G12V   binding  destab  <- basis: "WT-like RAF1 binding; impaired
                                       GTP hydrolysis"

All three were `neutral` in both the backlog and the canonical. The v7.6 review
(decision units D09 and D18, Gate B proposed, Gate C approved as
"source-cleared", confidence HIGH) changed them to `destab`. The bases were not
changed and still record intact binding, so the change was not supported by the
source it claimed to be cleared against.

For KRAS G12D and G12V the biochemistry is not in question: these are the
canonical GTPase-deficient oncogenic substitutions. Impaired intrinsic and
GAP-stimulated hydrolysis locks KRAS in the GTP-bound state and RAF1 engagement
is preserved or enhanced; the curated basis itself records that G12 sits 13.6 A
from the RAF1 interface. A binding-destabilizing expectation there is wrong.

For MSH2 C697F the basis records a catalytic ATPase defect with the protein
stable and the MSH6 interaction intact. The same reasoning applies: the lesion
is catalytic, which is not one of the three modeled axes.

This script reverses those three cells only. It does not touch the other 15
v7.6 decision units, which remain in force.

WHAT IT DOES NOT DO
-------------------
It does not make the benchmark agree with COMAVI. Two of the three reversals
move a variant into a class COMAVI grades correctly (G12D, G12V become
structurally_silent and grade `consistent`); the third moves MSH2 C697F into a
class COMAVI grades `inconsistent`, because COMAVI predicts 13.3 kcal/mol of
monomer destabilization on a buried-core residue the literature reports as
stable. That is a COMAVI false positive and the correction exposes it.

THE REPRODUCTION GATE
---------------------
Inherited from apply_ledger_corrections_v76.py: `scored_61var_canonical.csv` is
a frozen artifact that the v7.2-v7.5 stages mutated in place, so blind
re-derivation of Stage 6/6b does not reproduce it on four curated rows
(brca1_brct R1699L and R1699Q, vhl_elonginc W117R, troponin_ic R162W). This
script refuses to write any cell on a row re-derivation cannot reproduce. The
three corrections are disjoint from those four, so the gate passes.

Usage
-----
    python scripts/apply_ledger_corrections_v79.py [--dry-run] [--check]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

# Import the derivation path and cell renderer from v7.6 rather than copying
# them, so this correction cannot drift from the one that preceded it.
from apply_ledger_corrections_v76 import _cell, rederive  # noqa: E402

CANONICAL = ROOT / "reference_outputs" / "scored_61var_canonical.csv"
LEDGER = ROOT / "reference_outputs" / "COMAVI_evidence_ledger.csv"
SNAPDIR = ROOT / "reference_outputs"
MANIFEST = SNAPDIR / "COMAVI_v79_change_manifest.csv"
KEY = "variant"

# system, variant, column, from, to, the v7.6 decision unit being reversed
CORRECTIONS = [
    ("msh2_msh6", "C697F", "expected_ddg_binding", "destab", "neutral", "D18"),
    ("kras_craf", "G12D", "expected_ddg_binding", "destab", "neutral", "D09"),
    ("kras_craf", "G12V", "expected_ddg_binding", "destab", "neutral", "D09"),
]

REASON = ("v7.9: reversed -- the v7.6 token contradicted this row's own "
          "evidence_basis, which records intact binding")


def patch_ledger(dry_run: bool) -> int:
    """Bring the evidence ledger back into agreement with its own bases."""
    led = pd.read_csv(LEDGER)
    if "correction_note" not in led.columns:
        led["correction_note"] = ""
    n = 0
    for system, variant, column, frm, to, unit in CORRECTIONS:
        axis = column.replace("expected_ddg_", "")
        m = ((led["system"] == system) & (led["variant"] == variant)
             & (led["axis"] == axis))
        if m.sum() != 1:
            raise SystemExit("ABORT: ledger key %s/%s/%s matched %d rows"
                             % (system, variant, axis, m.sum()))
        stored = str(led.loc[m, "expected_token"].iloc[0])
        if stored != frm:
            raise SystemExit("ABORT: %s/%s/%s expected_token is %r, expected %r"
                             % (system, variant, axis, stored, frm))
        led.loc[m, "expected_token"] = to
        led.loc[m, "correction_note"] = "%s (%s)" % (REASON, unit)
        n += 1
    print("  ledger: %d tokens reverted to neutral" % n)
    if dry_run:
        return n
    ts = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    shutil.copy2(LEDGER, SNAPDIR / ("COMAVI_evidence_ledger.pre_v79_%s.csv" % ts))
    led.to_csv(LEDGER, index=False)
    print("  wrote    -> %s (%d rows)" % (LEDGER.name, len(led)))
    return n


def check() -> int:
    """Verify the correction is already applied and internally consistent."""
    base = pd.read_csv(CANONICAL, low_memory=False)
    led = pd.read_csv(LEDGER)
    bad = []
    for system, variant, column, frm, to, _unit in CORRECTIONS:
        axis = column.replace("expected_ddg_", "")
        row = base[(base["system"] == system) & (base[KEY] == variant)]
        if len(row) != 1:
            bad.append("%s/%s matched %d canonical rows" % (system, variant, len(row)))
            continue
        got = _cell(row[column].iloc[0])
        if got != to:
            bad.append("canonical %s/%s %s is %r, expected %r"
                       % (system, variant, column, got, to))
        lr = led[(led["system"] == system) & (led["variant"] == variant)
                 & (led["axis"] == axis)]
        if len(lr) == 1 and str(lr["expected_token"].iloc[0]) != to:
            bad.append("ledger %s/%s/%s is %r, expected %r"
                       % (system, variant, axis, lr["expected_token"].iloc[0], to))
    if bad:
        print("FAIL: v7.9 correction not applied or inconsistent")
        for b in bad:
            print("  - %s" % b)
        return 1
    print("PASS: all %d v7.9 corrections present in canonical and ledger"
          % len(CORRECTIONS))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--canonical", default=str(CANONICAL))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true",
                    help="verify the correction is applied; write nothing")
    args = ap.parse_args()

    if args.check:
        return check()

    path = Path(args.canonical)
    base = pd.read_csv(path, low_memory=False)
    print("v7.9: %d token reversals over %d variants"
          % (len(CORRECTIONS), len({c[1] for c in CORRECTIONS})))

    patched = base.copy()
    for system, variant, column, frm, to, _unit in CORRECTIONS:
        msk = (patched["system"] == system) & (patched[KEY] == variant)
        if msk.sum() != 1:
            raise SystemExit("ABORT: %s/%s matched %d rows" % (system, variant, msk.sum()))
        stored = _cell(patched.loc[msk, column].iloc[0])
        if stored != frm:
            raise SystemExit("ABORT: %s/%s %s is %r, expected %r -- refusing to write"
                             % (system, variant, column, stored, frm))
        patched.loc[msk, column] = to

    A = rederive(base)
    B = rederive(patched)
    shared = [c for c in A.columns if c in base.columns]

    moved = {c: A.index[A[c].map(_cell).values != B[c].map(_cell).values] for c in shared}
    moved = {c: ix for c, ix in moved.items() if len(ix)}
    moved_rows = sorted({base.at[i, KEY] for ix in moved.values() for i in ix})

    unreproducible = sorted({
        base.at[i, KEY]
        for c in shared
        for i in base.index[base[c].map(_cell).values != A[c].map(_cell).values]
    })
    collisions = sorted(set(moved_rows) & set(unreproducible))
    print("  re-derivation reproduces %d/%d variants; curated/unreproducible: %d %s"
          % (len(base) - len(unreproducible), len(base), len(unreproducible),
             unreproducible))
    if collisions:
        raise SystemExit("ABORT: correction would write to curated rows that "
                         "re-derivation cannot reproduce: %s" % collisions)
    print("  GATE PASS: %d cells across %d columns / %d variants; no curated-row collision"
          % (sum(len(ix) for ix in moved.values()), len(moved), len(moved_rows)))

    out = base.copy()
    manifest = []
    for system, variant, column, frm, to, unit in CORRECTIONS:
        msk = (out["system"] == system) & (out[KEY] == variant)
        manifest.append({"variant": variant, "system": system, "column": column,
                         "stored": frm, "written": to, "kind": "token",
                         "reverses_v76_unit": unit})
        out.loc[msk, column] = to
    for c, ix in moved.items():
        for i in ix:
            manifest.append({"variant": base.at[i, KEY], "system": base.at[i, "system"],
                             "column": c, "stored": _cell(base.at[i, c]),
                             "written": _cell(B.at[i, c]), "kind": "derived",
                             "reverses_v76_unit": ""})
            out.at[i, c] = B.at[i, c]

    man = pd.DataFrame(manifest)
    print("\n  what moves:")
    for v in sorted({m["variant"] for m in manifest}):
        sub = man[man.variant == v]
        cls = sub[sub.column == "expected_mech_class"]
        gr = sub[sub.column == "mech_consistency_t25"]
        print("    %-7s class %s -> %s   grade %s -> %s"
              % (v,
                 cls.stored.iloc[0] if len(cls) else "(unchanged)",
                 cls.written.iloc[0] if len(cls) else "",
                 gr.stored.iloc[0] if len(gr) else "(unchanged)",
                 gr.written.iloc[0] if len(gr) else ""))

    if args.dry_run:
        print("\n  --dry-run: nothing written")
        print(man.groupby("kind").size().to_string())
        patch_ledger(dry_run=True)
        return 0

    ts = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    snap = SNAPDIR / ("scored_61var_canonical.pre_v79_%s.csv" % ts)
    shutil.copy2(path, snap)
    out.to_csv(path, index=False)
    man.to_csv(MANIFEST, index=False)
    print("\n  snapshot -> %s" % snap.name)
    print("  wrote    -> %s  (%d rows, %d cols)" % (path.name, len(out), len(out.columns)))
    print("  manifest -> %s  (%d cells)" % (MANIFEST.name, len(man)))
    patch_ledger(dry_run=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
