#!/usr/bin/env python
"""Does COMAVI's structural agreement depend on how well-grounded the ground truth is?

A reviewer handed the evidence ledger will immediately ask whether the
headline agreement is propped up by weakly-grounded rows. This script answers
it, and the answer has a twist worth stating carefully.

CRUDE COMPARISON: agreement is LOWER on axes backed by quantitative
measurement than on axes backed by qualitative or inferred ground truth.
Read naively that says the pipeline does worse where the truth is harder --
which would be an alarming result.

IT IS CONFOUNDED. Quantitative measurements were sought precisely for the
destabilizers (the informative, hard-to-call cases); qualitative and inferred
tokens are overwhelmingly `neutral` (the easy call, and the one a
conservative predictor gets right by default). Evidence grade and expected
token are therefore strongly associated, and the crude contrast is a Simpson's
paradox artifact. Stratified by expected token, the association disappears.

This script emits both, plus the axis-availability accounting that explains
why the committed-axis count and the gradeable-axis count differ.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
LEDGER = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"
OUT_CSV = REPO / "reference_outputs" / "COMAVI_evidence_stratified_agreement.csv"
OUT_JSON = REPO / "reference_outputs" / "COMAVI_evidence_stratified_agreement.json"

T = 2.5  # canonical calling threshold
E1 = "E1_quantitative_energetic"
E2 = "E2_quantitative_functional"
QUANT = (E1, E2)
# structural_agreement_by_axis names the complex-fold axis "fold"; the
# evidence ledger uses the ground-truth column name "fold_complex".
AXIS_ALIAS = {"fold": "fold_complex"}


def mantel_haenszel(strata):
    """strata: list of (a, b, c, d) 2x2 counts. Returns (OR_MH, chi2, p)."""
    num = den = a_sum = e_sum = v_sum = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
        a_sum += a
        e_sum += (a + b) * (a + c) / n
        if n > 1:
            v_sum += ((a + b) * (c + d) * (a + c) * (b + d)) / (n ** 2 * (n - 1))
    if den == 0 or v_sum == 0:
        return float("nan"), float("nan"), float("nan")
    x2 = (abs(a_sum - e_sum) - 0.5) ** 2 / v_sum
    return num / den, x2, float(1 - chi2.cdf(x2, 1))


def rate_table(frame, by):
    t = frame.groupby(by).agg(ok=("n_ok", "sum"), n=("n_grade", "sum"))
    t["rate"] = (t.ok / t.n).round(3)
    return t


def main() -> int:
    df = pd.read_csv(CANON)
    led = pd.read_csv(LEDGER)
    partners = ac.discover_partners(df)

    rows = []
    for _, r in df.iterrows():
        per = ac.structural_agreement_by_axis(r, partners, T, T, T)
        for axis, (n_ok, n_grade) in per.items():
            rows.append(dict(system=r["system"], variant=r["variant"], axis=axis,
                             n_ok=n_ok, n_grade=n_grade))
    sa = pd.DataFrame(rows)
    sa["axis_led"] = sa.axis.replace(AXIS_ALIAS)

    merged = sa.merge(
        led[["system", "variant", "axis", "evidence_type",
             "evidence_directness", "expected_token"]],
        left_on=["system", "variant", "axis_led"],
        right_on=["system", "variant", "axis"],
        how="left", suffixes=("", "_led"),
    )
    grade = merged[merged.n_grade > 0].copy()
    ddg = grade[grade.axis != "tier"].copy()

    unmatched = int(ddg.evidence_type.isna().sum())
    if unmatched:
        print(f"{unmatched} gradeable ddG axes have no evidence-ledger row",
              file=sys.stderr)
        return 1

    ddg["quant"] = ddg.evidence_type.isin(QUANT)

    # ---- accounting: committed vs gradeable -------------------------------
    n_committed = len(led)
    n_ddg_gradeable = int(ddg.n_grade.sum())
    n_gated_out = n_committed - n_ddg_gradeable

    # ---- crude contrast ---------------------------------------------------
    crude = rate_table(ddg, "quant")
    a = int(crude.loc[True, "ok"]); b = int(crude.loc[True, "n"] - a)
    c = int(crude.loc[False, "ok"]); d = int(crude.loc[False, "n"] - c)
    or_crude, p_crude = fisher_exact([[a, b], [c, d]])

    # ---- confound: token composition differs by evidence grade -----------
    comp = pd.crosstab(ddg.quant, ddg.expected_token,
                       values=ddg.n_grade, aggfunc="sum").fillna(0).astype(int)

    # ---- stratified contrast ---------------------------------------------
    strata, per_tok = [], {}
    for tok in sorted(ddg.expected_token.dropna().unique()):
        t = rate_table(ddg[ddg.expected_token == tok], "quant")
        if t.shape[0] < 2:
            per_tok[tok] = {"note": "single evidence-grade stratum",
                            **{str(k): int(v) for k, v in
                               t.reset_index().iloc[0][["ok", "n"]].items()}}
            continue
        aa = int(t.loc[True, "ok"]); bb = int(t.loc[True, "n"] - aa)
        cc = int(t.loc[False, "ok"]); dd = int(t.loc[False, "n"] - cc)
        strata.append((aa, bb, cc, dd))
        o, p = fisher_exact([[aa, bb], [cc, dd]])
        per_tok[tok] = dict(quant_ok=aa, quant_n=aa + bb,
                            soft_ok=cc, soft_n=cc + dd,
                            quant_rate=round(aa / (aa + bb), 3),
                            soft_rate=round(cc / (cc + dd), 3),
                            OR=round(float(o), 3), p=round(float(p), 4))
    or_mh, x2_mh, p_mh = mantel_haenszel(strata)

    by_type = rate_table(ddg, "evidence_type")
    by_dir = rate_table(ddg, "evidence_directness")
    ddg.to_csv(OUT_CSV, index=False)

    summary = dict(
        threshold=T,
        accounting=dict(
            n_committed_axes=n_committed,
            n_ddg_axes_gradeable=n_ddg_gradeable,
            n_committed_but_gated_out=n_gated_out,
            n_tier_axes_gradeable=int(grade[grade.axis == "tier"].n_grade.sum()),
            n_total_gradeable=int(grade.n_grade.sum()),
            note=("Committed axes carry ground truth; a committed axis is "
                  "gradeable only if the prediction also clears confidence "
                  "and internal-CI gating. The difference is the cost of "
                  "those gates, not missing curation."),
        ),
        crude=dict(
            quant_ok=a, quant_n=a + b, quant_rate=round(a / (a + b), 3),
            soft_ok=c, soft_n=c + d, soft_rate=round(c / (c + d), 3),
            OR=round(float(or_crude), 3), p=round(float(p_crude), 4),
        ),
        token_composition_by_grade={
            ("quantitative" if k else "soft"): v
            for k, v in comp.to_dict("index").items()
        },
        stratified=dict(
            per_token=per_tok,
            OR_mantel_haenszel=round(float(or_mh), 3),
            chi2=round(float(x2_mh), 3),
            p=round(float(p_mh), 4),
        ),
        by_evidence_type={k: dict(ok=int(v["ok"]), n=int(v["n"]),
                                  rate=float(v["rate"]))
                          for k, v in by_type.to_dict("index").items()},
        by_directness={k: dict(ok=int(v["ok"]), n=int(v["n"]),
                               rate=float(v["rate"]))
                       for k, v in by_dir.to_dict("index").items()},
        interpretation=(
            "The crude contrast is confounded: quantitative measurement was "
            "sought for destabilizers, while qualitative and inferred tokens "
            "are predominantly neutral. Stratified by expected token the "
            "association is absent, so evidence grade does not explain the "
            "headline agreement."
        ),
    )
    OUT_JSON.write_text(json.dumps(summary, indent=1))

    print(f"committed axes {n_committed} | ddG gradeable {n_ddg_gradeable} "
          f"| gated out {n_gated_out} | tier gradeable "
          f"{summary['accounting']['n_tier_axes_gradeable']}")
    print(f"\ncrude: quant {a}/{a+b} = {a/(a+b):.3f} vs soft "
          f"{c}/{c+d} = {c/(c+d):.3f}  OR={or_crude:.2f} p={p_crude:.4f}")
    print("\ntoken composition by evidence grade (gradeable ddG axes):")
    print(comp.to_string())
    print(f"\nstratified by expected token: OR_MH={or_mh:.2f} "
          f"chi2={x2_mh:.2f} p={p_mh:.3f}")
    print("\nby evidence type:")
    print(by_type.to_string())
    print(f"\nwrote {OUT_CSV.relative_to(REPO)}, {OUT_JSON.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
