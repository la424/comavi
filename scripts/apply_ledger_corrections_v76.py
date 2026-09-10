#!/usr/bin/env python
"""
v7.6 -- Apply the approved evidence-ledger ground-truth corrections.

In-place mutator on the canonical scoring table, following the same
convention as apply_canonical_corrections_v72 / apply_brct_pooling_v73 /
apply_brct_annotations_v74 / apply_brct_truth_threshold_v75: snapshot
first, mutate in place, print a manifest.

What it changes
---------------
The `expected_*` ground-truth tokens for the decision units approved in the
Gate B / Gate C / Gate D review, then re-derives the columns downstream of
those tokens (Stage 6 / 6b of apply_concordance_v5.py).

The reproduction gate
---------------------
`scored_61var_canonical.csv` is a frozen artifact: the v7.2-v7.5 stages
mutated it in place after apply_concordance_v5.py last ran, so a blind
re-derivation of Stage 6/6b does NOT reproduce it everywhere. Re-derivation
reproduces 47 of 61 variants exactly; it does not reproduce the 12
brca1_brct rows (curated by v7.3/v7.4/v7.5) nor W117R (vhl_elonginc) and
R162W (troponin_ic) (curated by v7.2).

This script therefore refuses to write any cell on a row that
re-derivation cannot reproduce. The approved corrections happen to be
disjoint from those 14 rows, so the gate passes; if a future correction
touches one of them, this script will abort rather than silently revert a
curated value.

Usage
-----
    python scripts/apply_ledger_corrections_v76.py [--dry-run]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import apply_concordance_v5 as ac  # noqa: E402

CANONICAL = ROOT / "reference_outputs" / "scored_61var_canonical.csv"
SCENARIO = ROOT / "inputs" / "COMAVI_phase5_scenario_definition.csv"
SNAPDIR = ROOT / "reference_outputs"
KEY = "variant"

_NA = {"", "nan", "none", "n/a", "na", "<na>", "null"}


def _cell(v) -> str:
    """Dtype-insensitive scalar rendering, so float 1.0 == int 1 == '1'."""
    if v is None:
        return ""
    if isinstance(v, float) and np.isnan(v):
        return ""
    if isinstance(v, (int, float, np.integer, np.floating)):
        return f"{float(v):.9g}"
    s = str(v).strip()
    return "" if s.lower() in _NA else s


def rederive(m: pd.DataFrame) -> pd.DataFrame:
    """Replicate Stage 6 + 6b of apply_concordance_v5.py (the token-dependent block)."""
    m = m.copy()
    P = ac.discover_partners(m)
    m["expected_mech_class"] = m.apply(ac.derive_expected_mech_class, axis=1)
    axes_by_idx = {i: ac.classify_axis_status(r) for i, r in m.iterrows()}

    def _grade_one(r, mech_col):
        g, fp, mp = ac.grade_mechanism_consistency(
            r, r.get(mech_col), r.get("expected_mech_class"), axes_by_idx.get(r.name)
        )
        return g, ",".join(fp) if fp else "", ",".join(mp) if mp else ""

    tags = [t for t, _ in ac.THRESHOLD_SPECS]
    for suf in tags:
        for pre, col in (("", f"mech_{suf}"), ("nbhd_", f"nbhd_mech_{suf}")):
            g = [_grade_one(r, col) for _, r in m.iterrows()]
            m[f"{pre}mech_consistency_{suf}"] = [x[0] for x in g]
            m[f"{pre}mech_false_positive_axes_{suf}"] = [x[1] for x in g]
            m[f"{pre}mech_missed_positive_axes_{suf}"] = [x[2] for x in g]

    for pre in ("", "nbhd_"):
        m[f"{pre}mech_consistency_summary"] = m[f"{pre}mech_consistency_t25"]
        m[f"{pre}mech_consistency_threshold_stable"] = m.apply(
            lambda r, _p=pre: len({r.get(f"{_p}mech_consistency_{s}") for s in tags}) == 1,
            axis=1,
        )

    for tag, t in ac.THRESHOLD_SPECS:
        mt, ft, bt = (
            (t["monomer"], t["fold"], t["binding"]) if isinstance(t, dict) else (float(t),) * 3
        )
        sa = m.apply(lambda r: ac.compute_structural_agreement(r, P, mt, ft, bt), axis=1)
        m[f"structural_agreement_n_{tag}"] = [r[0] for r in sa]
        m[f"structural_agreement_d_{tag}"] = [r[1] for r in sa]
        m[f"structural_agreement_{tag}"] = [r[0] / r[1] if r[1] > 0 else None for r in sa]

        da = m.apply(lambda r: ac.compute_directional_agreement(r, P, mt, ft, bt), axis=1)
        m[f"directional_agreement_full_{tag}"] = [r[0] for r in da]
        m[f"directional_agreement_half_{tag}"] = [r[1] for r in da]
        m[f"directional_agreement_d_{tag}"] = [r[2] for r in da]
        m[f"directional_agreement_{tag}"] = [
            (r[0] + 0.5 * r[1]) / r[2] if r[2] > 0 else None for r in da
        ]

    vr = m.apply(lambda r: ac.compute_axis_votes(r, P), axis=1)
    for c in (
        "ddg_monomer_vote_strict", "ddg_monomer_vote_relaxed",
        "ddg_fold_vote_strict", "ddg_fold_vote_relaxed",
        "ddg_binding_vote_strict", "ddg_binding_vote_relaxed",
    ):
        m[c] = vr.apply(lambda d, _c=c: d.get(_c))

    m["tier_structural_signal_type"] = m.apply(
        lambda r: ac.compute_tier_structural_signal_type(r, P), axis=1
    )
    m["evaluation_note"] = m.apply(
        lambda r: ac.compute_evaluation_note(
            r,
            (r.get("structural_agreement_n_t25", 0) or 0,
             r.get("structural_agreement_d_t25", 0) or 0),
            r.get("mech_consistency_t25"),
        ),
        axis=1,
    )
    return m


LEDGER = ROOT / "reference_outputs" / "COMAVI_evidence_ledger.csv"
WITHDRAWN = ROOT / "reference_outputs" / "COMAVI_evidence_ledger_withdrawn_v76.csv"


def patch_ledger(changes: pd.DataFrame, dry_run: bool) -> None:
    """Bring the evidence ledger into agreement with the corrected canonical.

    The ledger enumerates *committed* expected tokens. A correction to
    `unknown` withdraws the commitment, so the row leaves the ledger; its
    evidence fields are preserved in COMAVI_evidence_ledger_withdrawn_v76.csv
    rather than discarded, so the audit trail survives the withdrawal.
    """
    led = pd.read_csv(LEDGER)
    chg = changes.assign(axis=changes["column"].str.replace("expected_ddg_", "", regex=False))
    keys = led.set_index(["system", "variant", "axis"]).index

    withdraw = chg[chg["final_token"].fillna("unknown") == "unknown"]
    retoken = chg[chg["final_token"].fillna("unknown") != "unknown"]

    wmask = keys.isin(list(withdraw.set_index(["system", "variant", "axis"]).index))
    wrows = led[wmask].copy()
    wrows["withdrawn_by"] = "v7.6"
    wrows["withdrawal_reason"] = withdraw.set_index(["system", "variant", "axis"]).reindex(
        pd.MultiIndex.from_frame(wrows[["system", "variant", "axis"]])
    )["final_source"].values

    out = led[~wmask].copy()
    n_ret = 0
    for _, r in retoken.iterrows():
        m = ((out["system"] == r["system"]) & (out["variant"] == r["variant"])
             & (out["axis"] == r["axis"]))
        if m.sum() != 1:
            raise SystemExit(f"ABORT: ledger key {r['system']}/{r['variant']}/{r['axis']} "
                             f"matched {m.sum()} rows")
        out.loc[m, "expected_token"] = r["final_token"]
        n_ret += 1

    print(f"  ledger: {len(led)} rows -> {len(out)} "
          f"({len(wrows)} commitments withdrawn, {n_ret} retokened)")
    if dry_run:
        return
    ts = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    shutil.copy2(LEDGER, SNAPDIR / f"COMAVI_evidence_ledger.pre_v76_{ts}.csv")
    out.to_csv(LEDGER, index=False)
    wrows.to_csv(WITHDRAWN, index=False)
    print(f"  wrote    -> {LEDGER.name} ({len(out)} rows)")
    print(f"  withdrawn-> {WITHDRAWN.name} ({len(wrows)} rows, evidence preserved)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--canonical", default=str(CANONICAL))
    ap.add_argument("--scenario", default=str(SCENARIO))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--ledger-only", action="store_true",
                    help="patch only the evidence ledger (canonical already corrected)")
    args = ap.parse_args()

    if args.ledger_only:
        s = pd.read_csv(args.scenario)
        patch_ledger(s[s["is_token_change"].astype(bool)], args.dry_run)
        return 0

    path = Path(args.canonical)
    base = pd.read_csv(path, low_memory=False)
    scen = pd.read_csv(args.scenario)
    changes = scen[scen["is_token_change"].astype(bool)]
    print(f"v7.6: {len(changes)} approved token changes over "
          f"{changes['variant'].nunique()} variants, {changes['system'].nunique()} systems")

    patched = base.copy()
    for _, r in changes.iterrows():
        msk = (patched["system"] == r["system"]) & (patched[KEY] == r["variant"])
        if msk.sum() != 1:
            raise SystemExit(f"ABORT: {r['system']}/{r['variant']} matched {msk.sum()} rows")
        patched.loc[msk, r["column"]] = r["final_token"]

    A = rederive(base)      # re-derived from unpatched tokens
    B = rederive(patched)   # re-derived from corrected tokens
    shared = [c for c in A.columns if c in base.columns]

    # cells the correction actually moves
    moved = {c: A.index[A[c].map(_cell).values != B[c].map(_cell).values] for c in shared}
    moved = {c: ix for c, ix in moved.items() if len(ix)}
    moved_rows = sorted({base.at[i, KEY] for ix in moved.values() for i in ix})

    # gate: re-derivation must reproduce the frozen table on every row we touch
    unreproducible = sorted({
        base.at[i, KEY]
        for c in shared
        for i in base.index[base[c].map(_cell).values != A[c].map(_cell).values]
    })
    collisions = sorted(set(moved_rows) & set(unreproducible))
    print(f"  re-derivation reproduces {len(base) - len(unreproducible)}/{len(base)} variants; "
          f"curated/unreproducible: {len(unreproducible)}")
    if collisions:
        raise SystemExit(
            "ABORT: correction would write to curated rows that re-derivation "
            f"cannot reproduce: {collisions}"
        )
    print(f"  GATE PASS: {sum(len(ix) for ix in moved.values())} cells to write across "
          f"{len(moved)} columns / {len(moved_rows)} variants; no collision with curated rows")

    out = base.copy()
    manifest = []
    for _, r in changes.iterrows():
        msk = (out["system"] == r["system"]) & (out[KEY] == r["variant"])
        manifest.append({"variant": r["variant"], "system": r["system"], "column": r["column"],
                         "stored": _cell(out.loc[msk, r["column"]].iloc[0]),
                         "written": r["final_token"], "kind": "token"})
        out.loc[msk, r["column"]] = r["final_token"]
    for c, ix in moved.items():
        for i in ix:
            manifest.append({"variant": base.at[i, KEY], "system": base.at[i, "system"],
                             "column": c, "stored": _cell(base.at[i, c]),
                             "written": _cell(B.at[i, c]), "kind": "derived"})
            out.at[i, c] = B.at[i, c]

    man = pd.DataFrame(manifest)
    if args.dry_run:
        print("  --dry-run: nothing written")
        print(man.groupby("kind").size().to_string())
        return 0

    ts = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    snap = SNAPDIR / f"scored_61var_canonical.pre_v76_{ts}.csv"
    shutil.copy2(path, snap)
    out.to_csv(path, index=False)
    man_path = SNAPDIR / "COMAVI_v76_change_manifest.csv"
    man.to_csv(man_path, index=False)
    print(f"  snapshot -> {snap.name}")
    print(f"  wrote    -> {path.name}  ({len(out)} rows, {len(out.columns)} cols)")
    print(f"  manifest -> {man_path.name}  ({len(man)} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
