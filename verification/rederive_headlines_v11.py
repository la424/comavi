#!/usr/bin/env python3
"""Re-derive the manuscript's headline values from the canonical table alone.

WHY THIS EXISTS
---------------
Every gate in this repository compares a record with the generator that wrote
it, or prose with a record. None of them would notice if a generator and its
record were wrong together. This script is written from the published
definitions -- Table 3 (the grading rubric), Table 2 (the metrics) and the
Methods (confidence gate, priority-score formula, direction agreement) -- and
imports no generator. It reads only primitive inputs: the canonical table
(energies, model confidence, tiers, stored calls), the evidence ledger
(committed expectations and the assays behind them), the measured-energy
comparisons and the labelled pathogenicity set. Each re-derived value is then
looked for, in the manuscript's own phrasing, in the extracted manuscript text.

The chain is checked link by link:
  energies -> calls      the stored mechanism calls follow from the energies under
                         the Methods gate (site pLDDT >= 50) and the threshold;
  calls -> grades        every stored grade follows from Table 3;
  grades -> headlines    scores, arms, rates, structure-blind baselines, margin,
                         direction agreement, resampling and permutation results;
  energies -> priority   the priority score follows from its published formula.

Stochastic quantities (system resampling, label permutation) are re-derived
with independent seeds and draw counts and compared within a stated Monte
Carlo tolerance; they are not expected to reproduce digit for digit.

Two findings of the first run are encoded rather than hidden:
  * Direction agreement admits an energy axis only when the interval of its
    prediction (replicate mean +/- 1.96 replicate SD) excludes zero. The
    published counts reproduce only with this rule, which the Methods did not
    state until v11; all 18 of the 97 expectations it removes are removed by
    this check, none by model confidence.
  * Measured-effect recovery at the prediction-interval limit used 2.9 kcal/mol
    for binding comparisons instead of 3.5, so the manuscript said 11 of 42
    where the per-axis definition gives 8 (fixed in build_measured_effect_tables).
  * The label-permutation null and the variant bootstrap ran over all 61 scored
    rows, mixing the calls of excluded candidates into the null; both now run on
    the 56 benchmark variants (fixed in verification/stress_tests.py).
  * The call vocabulary has no label for an isolated subunit stabilized while
    the complex is destabilized, so PIK3CA H1047R (isolated subunit -1.59
    kcal/mol) is called on the complex axis alone at 1.0 and 1.5 kcal/mol.
    It is a non-disruptive control graded inconsistent either way, and the
    isolated-subunit rows are computed from energies, so no reported value
    depends on it. Any other call that does not follow from the energies fails.

Usage:
    python verification/rederive_headlines_v11.py [--report path.csv]
"""
from __future__ import annotations

import argparse
import collections
import json
import pathlib
import re
import sys
from fractions import Fraction
from math import comb

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"
sys.path.insert(0, str(REPO / "scripts"))
import benchmark_population as bp  # noqa: E402  (the curated exclusion list only)

TAGS = ("t10", "t15", "t20", "t25", "tSAP")
THR = {"t10": (1.0, 1.0, 1.0), "t15": (1.5, 1.5, 1.5), "t20": (2.0, 2.0, 2.0), "t25": (2.5, 2.5, 2.5),
       "tSAP": (2.9, 2.9, 3.5)}                    # isolated subunit, assembled complex, binding
AXES = ("monomer", "fold_complex", "binding")
TAU = {"monomer": 2.9, "fold_complex": 2.9, "binding": 3.5}  # priority-score reference values
GM = {"consistent": Fraction(1), "partial": Fraction(1, 2), "inconsistent": Fraction(0)}
SILENT = "structurally_silent"
KNOWN_CALL_EXCEPTIONS = {("pi3k", "H1047R", "t10"), ("pi3k", "H1047R", "t15")}


# ------------------------------------------------------------------ labels
def label_axes(label):
    """What a call label asserts on each axis, parsed from its name."""
    if label.startswith(("No structural effect", "Interface variant", "Structural variant")):
        return {}
    if label == "PPI conflicting (mixed partner signals)":
        # partners move in opposite directions; some partner moves in each
        # direction, so a binding change of either sign counts as recovered.
        # It never occurs among the benchmark's calls, only as a constant.
        return {"binding": "mixed"}
    out = {}
    fold_axes = {"Monomer": ("monomer",), "Multimer": ("fold_complex",), "Both": ("monomer", "fold_complex")}
    ppi = re.search(r"PPI (de)?stabilization", label)
    if ppi:
        out["binding"] = "destab" if ppi.group(1) else "stab"
    m = re.match(r"(Monomer|Multimer|Both) fold (destab|stab)\w*\.?", label)
    shared = re.match(r"(Monomer|Multimer|Both) fold \+ PPI", label)   # one direction for fold and binding
    if m:
        for ax in fold_axes[m.group(1)]:
            out[ax] = "destab" if m.group(2) == "destab" else "stab"
    elif shared:
        for ax in fold_axes[shared.group(1)]:
            out[ax] = out["binding"]
    if not out:
        raise ValueError("unparsed call label: %r" % label)
    return out


# ------------------------------------------------------------------ inputs
CAN = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
B = bp.benchmark(CAN).reset_index(drop=True)
GRADED = CAN[CAN.mech_consistency_t25.notna()]
assert set(zip(B.system, B.variant)) == set(zip(GRADED.system, GRADED.variant)), "benchmark != graded rows"
MECH = (B.expected_mech_class != SILENT).values
LED = pd.read_csv(RO / "COMAVI_evidence_ledger.csv", keep_default_na=False, dtype=str)
EXP = collections.defaultdict(dict)
ASSAY = {}
for r in LED.itertuples():
    EXP[(r.system, r.variant)][r.axis] = r.expected_token
    if r.complex_amount_assay:
        ASSAY[(r.system, r.variant)] = r.complex_amount_assay


def either_scope(key):
    """The evidence-resolution rule, from the ledger: a binding change measured only
    as an amount of complex, with assembled-complex stability untested."""
    e = EXP.get(key, {})
    return key in ASSAY and e.get("binding") in ("destab", "stab") and "fold_complex" not in e


EITHER = [either_scope(k) for k in zip(B.system, B.variant)]
stored_scope = (B.assembled_expectation_scope.fillna("") == "either_axis").tolist()
assert EITHER == stored_scope, "evidence-resolution scope does not follow from the ledger"
PARTNERS = sorted({m.group(1) for c in CAN.columns for m in [re.fullmatch(r"ddg_binding_([a-z0-9]+(?:_\d)?)", c)] if m})


def plddt_ok(row, ax, p):
    v = row.monomer_plddt if ax == "monomer" else row.get("multi_%s_plddt" % re.sub(r"_\d$", "", p))
    return pd.isna(v) or float(v) >= 50  # experimental structures carry no pLDDT


def interval_ok(row, ax, p):
    col = ("ddg_monomer_distinguishable_internal_from_0" if ax == "monomer" else
           "ddg_%s_%s_distinguishable_internal_from_0" % ("fold" if ax == "fold_complex" else "binding", p))
    v = row.get(col)
    return bool(v) if pd.notna(v) else False


def values(row, ax, gate):
    if ax == "monomer":
        return [float(row.ddg_monomer)] if pd.notna(row.ddg_monomer) and gate(row, "monomer", None) else []
    pre = "fold" if ax == "fold_complex" else "binding"
    return [float(row["ddg_%s_%s" % (pre, p)]) for p in PARTNERS
            if pd.notna(row.get("ddg_%s_%s" % (pre, p))) and gate(row, ax, p)]


# ------------------------------------------------------------------ energies -> calls
def calls_from_energy(row, t):
    th = dict(zip(AXES, THR[t]))
    out = {}
    for ax in AXES:
        v = values(row, ax, plddt_ok)
        if any(x >= th[ax] for x in v):
            out[ax] = "destab"
        elif any(x <= -th[ax] for x in v):
            out[ax] = "stab"
    return out


# ------------------------------------------------------------------ calls -> grades (Table 3)
def grade(fired, exp, cls, either):
    if cls == SILENT:
        return "consistent" if not fired else "inconsistent"
    change = {a: d for a, d in exp.items() if d in ("destab", "stab")}
    neutral = {a for a, d in exp.items() if d == "neutral"}
    recovered, opposite = set(), False
    for a, d in change.items():
        cands = [a] + (["fold_complex"] if either and a == "binding" else [])
        dirs = [fired[c] for c in cands if c in fired]
        if d in dirs or "mixed" in dirs:
            recovered.add(a)
        elif dirs:
            opposite = True
    if opposite or not recovered:
        return "inconsistent"
    unsupported = {a for a in fired if a in neutral}   # untested axes are neither credited nor penalized
    missed = set(change) - recovered
    if missed:
        return "inconsistent" if unsupported else "partial"
    return "partial" if unsupported else "consistent"


def score(labels_by_row, profiles):
    return sum(GM[grade(label_axes(L), *prof)] for L, prof in zip(labels_by_row, profiles))


PROFILES = [(EXP.get(k, {}), c, e) for k, c, e in zip(zip(B.system, B.variant), B.expected_mech_class, EITHER)]


# ------------------------------------------------------------------ direction agreement
def agreement(row, t, exp, either):
    th = dict(zip(AXES, THR[t]))
    gate = lambda r, a, p: plddt_ok(r, a, p) and interval_ok(r, a, p)  # noqa: E731

    def operative(ax):
        v = values(row, ax, gate)
        return max(v, key=abs) if v else None

    def ok(v, thr, d):
        fired = abs(v) >= thr
        return (fired and v > 0) if d == "destab" else (fired and v < 0) if d == "stab" else not fired
    out = {}
    for a, d in exp.items():
        if either and a == "binding":
            sign = 1 if d == "destab" else -1
            vs = [(operative(x), th[x]) for x in ("fold_complex", "binding")]
            vs = [(v, tt) for v, tt in vs if v is not None]
            if vs:
                same = any(abs(v) >= tt and v * sign > 0 for v, tt in vs)
                opp = any(abs(v) >= tt and v * sign < 0 for v, tt in vs)
                out["either"] = (int(same and not opp), 1)
            continue
        v = operative(a)
        if v is not None:
            out[a] = (int(ok(v, th[a], d)), 1)
    if pd.notna(row.comavi_tier):
        out["tier"] = (int((row.comavi_tier in ("Tier 1", "Tier 2")) == (row.expected_mech_class != SILENT)), 1)
    return out


# ------------------------------------------------------------------ priority score
def priority(row):
    ratios = [abs(v) / TAU[ax] for ax in AXES for v in values(row, ax, plddt_ok)]
    R = max(ratios) if ratios else 0.0
    T = int(str(row.comavi_tier).split()[-1])
    return R / (1 + R), (4 - T) / 3


def exact_perm_p(system, truth, strong):
    pop = pd.DataFrame({"system": list(system), "truth": list(truth), "strong": list(strong)})
    obs = int((pop["truth"] & pop["strong"]).sum())
    tot = {0: Fraction(1)}
    for _, g in pop.groupby("system"):
        n, k, s = len(g), int(g["truth"].sum()), int(g["strong"].sum())
        dist = {x: Fraction(comb(s, x) * comb(n - s, k - x), comb(n, k)) for x in range(max(0, k - (n - s)), min(k, s) + 1)}
        new = {}
        for a, pa in tot.items():
            for x, px in dist.items():
                new[a + x] = new.get(a + x, 0) + pa * px
        tot = new
    return sum(p for v, p in tot.items() if v >= obs)


def r3(x, nd=3):
    """Round once, half to even, on the exact value."""
    q = x if isinstance(x, Fraction) else Fraction(repr(float(x))).limit_denominator(10 ** 6)
    neg, q = q < 0, abs(q)
    s = q * 10 ** nd
    lo = s.numerator // s.denominator
    if s - lo > Fraction(1, 2) or (s - lo == Fraction(1, 2) and lo % 2):
        lo += 1
    out = str(lo).rjust(nd + 1, "0")
    out = out[:-nd] + "." + out[-nd:] if nd else out
    return ("\u2212" if neg and lo else "") + out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--report", default=None, help="write every comparison to this CSV")
    ap.add_argument("--draws", type=int, default=20000, help="resampling and permutation draws")
    args = ap.parse_args()
    T = (REPO / "docs" / "COMAVI_manuscript_current.txt").read_text()
    body = re.sub(r"\s+", " ", T[:T.find("[SI S1_Text]")] if "[SI S1_Text]" in T else T)
    rows, fails = [], []

    def claim(name, value, phrase):
        found = phrase in body
        rows.append({"quantity": name, "re_derived": value, "manuscript_phrase": phrase, "found": found})
        if not found:
            fails.append("%s: re-derived %s, but the manuscript does not say %r" % (name, value, phrase))

    def check(name, ok, detail):
        rows.append({"quantity": name, "re_derived": detail, "manuscript_phrase": "", "found": bool(ok)})
        if not ok:
            fails.append("%s: %s" % (name, detail))

    # 1. energies -> calls
    bad_calls = [(r.system, r.variant, t) for _, r in B.iterrows() for t in TAGS
                 if calls_from_energy(r, t) != label_axes(r["mech_%s" % t])]
    check("calls follow from energies (Methods gate)", set(bad_calls) <= KNOWN_CALL_EXCEPTIONS,
          "%d of %d calls reproduce; exceptions %s" % (len(B) * 5 - len(bad_calls), len(B) * 5, sorted(bad_calls)))
    # 2. calls -> grades
    bad_grades = [(r.system, r.variant, t) for (_, r), prof in zip(B.iterrows(), PROFILES) for t in TAGS
                  if grade(label_axes(r["mech_%s" % t]), *prof) != r["mech_consistency_%s" % t]]
    check("grades follow from Table 3", not bad_grades, "%d of %d grades reproduce" % (len(B) * 5 - len(bad_grades), len(B) * 5))
    # 3. scores, arms and rates
    tot = {t: score(B["mech_%s" % t], PROFILES) for t in TAGS}
    s25 = tot["t25"]
    claim("mechanism-pattern score at 2.5", r3(s25 / len(B)), "the mechanism-pattern score was %s of 56 (%s)"
          % (r3(s25, 1), r3(s25 / len(B))))
    det = att = rej = 0
    for (_, r), (exp, cls, either) in zip(B.iterrows(), PROFILES):
        f = label_axes(r.mech_t25)
        if cls == SILENT:
            rej += not f
            continue
        impl = {a for a, d in exp.items() if d in ("destab", "stab")} | ({"fold_complex"} if either else set())
        det += bool(f)
        att += bool(set(f) & impl)
    n_m, n_s = int(MECH.sum()), int((~MECH).sum())
    claim("detection", det, "reached the threshold for %d of the %d mechanisms (detection, %s)" % (det, n_m, r3(Fraction(det, n_m))))
    claim("attribution", att, "reached it for %d (attribution, %s)" % (att, r3(Fraction(att, n_m))))
    claim("correct rejection", rej, "No axis reached it for %d of the %d non-disruptive variants (correct rejection, %s)"
          % (rej, n_s, r3(Fraction(rej, n_s))))
    # 3b. what the assembled axes add, by class (counts the colleague review asked for; v8.9 #25-#29)
    WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
    only = [(r.expected_mech_class, r.system) for _, r in B.iterrows()
            if r.expected_mech_class != SILENT and label_axes(r.mech_t25) and "monomer" not in label_axes(r.mech_t25)]
    by = collections.Counter(c for c, _ in only)
    tot_cls = collections.Counter(B.expected_mech_class[MECH])
    assert "brca1_brct" in set(B.system)          # the one single-chain system
    n_cx = B[B.system != "brca1_brct"].system.nunique()
    claim("mechanisms found only through the assembled axes, by class", "%d %s" % (len(only), dict(by)),
          "The %d mechanisms found only through the assembled axes are %d of the %d interface destabilizers, %d of the %d "
          "mixed mechanisms and %s of the %d fold mechanisms, drawn from %d of the %d complexes"
          % (len(only), by["ppi_destab_mechanism"], tot_cls["ppi_destab_mechanism"], by["mixed_structural"],
             tot_cls["mixed_structural"], WORDS.get(by["fold_mechanism"], str(by["fold_mechanism"])),
             tot_cls["fold_mechanism"], len({sy for _, sy in only}), n_cx))
    # 4. structure-blind strategies
    vocab = sorted({L for t in TAGS for L in CAN["mech_%s" % t].dropna()} |
                   set(json.loads((RO / "COMAVI_rubric_baselines.json").read_text())["every_constant_label"]["labels"]))
    vocab = [L for L in vocab if not L.startswith("Structure ")]
    const = {L: score([L] * len(B), PROFILES) for L in vocab}
    best = max(const.values())
    check("constant-label vocabulary", len(vocab) == 31, "%d graded constant labels" % len(vocab))
    claim("strongest structure-blind strategy", r3(best / len(B)), "against %s for the best structure-blind strategy" % r3(best / len(B)))
    claim("predict disruption on all three axes", r3(const["Both fold + PPI destabilization"] / len(B)),
          "all three axes scores %s" % r3(const["Both fold + PPI destabilization"] / len(B)))
    claim("margin over the strongest strategy", r3((s25 - best) / len(B)),
          "exceeds the best structure-blind strategy by %s" % r3((s25 - best) / len(B)))
    check("margin positive at every prespecified threshold",
          all(tot[t] > max(score([L] * len(B), PROFILES) for L in vocab) for t in TAGS), str({t: float(tot[t]) for t in TAGS}))
    # 4b. the Limitations paragraph: evidence composition, system concentration and the
    # dependence of localization on evidence resolution (added in v11 at review)
    WORDS_ = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
    classes = [{t_[:2] for t_ in LED[(LED.system == s_) & (LED.variant == v_)].evidence_type} for s_, v_ in zip(B.system, B.variant)]
    quant = [bool(c) and c <= {"E1", "E2"} for c in classes]          # every committed expectation quantitative
    qi = [i for i, q_ in enumerate(quant) if q_]
    g_of = lambda L, i: GM[grade(label_axes(L), *PROFILES[i])]        # noqa: E731
    q_comavi = sum(g_of(B.mech_t25[i], i) for i in qi) / len(qi)
    q_all = sum(g_of("Both fold + PPI destabilization", i) for i in qi) / len(qi)
    q_mech = [i for i in qi if MECH[i]]
    q_nd = [i for i in qi if not MECH[i]]
    arm = lambda L_of, idx: sum(g_of(L_of(i), i) for i in idx) / len(idx)  # noqa: E731
    q_bal = (arm(lambda i: B.mech_t25[i], q_mech) + arm(lambda i: B.mech_t25[i], q_nd)) / 2
    ceiling = max((arm(lambda i, L=L: L, q_mech) + arm(lambda i, L=L: L, q_nd)) / 2 for L in vocab)
    contrib = collections.defaultdict(Fraction)
    for i in range(len(B)):
        contrib[B.system[i]] += g_of(B.mech_t25[i], i) - g_of("No structural effect detected", i)
    att_no_rule = 0
    for (_, r), (exp, cls, either) in zip(B.iterrows(), PROFILES):
        if cls != SILENT:
            att_no_rule += bool(set(label_axes(r.mech_t25)) & {a for a, d in exp.items() if d in ("destab", "stab")})
    n_either = sum(1 for (exp, cls, either) in PROFILES if either and exp.get("binding") in ("destab", "stab"))
    claim("Limitations: quantitative subset, contributions and evidence resolution",
          "%d (%d non-disruptive), %s vs %s, balanced %s / %s; %s points, BRCT %s, Hb %s; %d vs %d"
          % (len(qi), len(q_nd), r3(q_all), r3(q_comavi), r3(q_bal), r3(ceiling, 1), r3(sum(contrib.values()), 1),
             r3(contrib["brca1_brct"], 1), r3(contrib["hemoglobin_tetramer"], 1), att_no_rule, att),
          "On the %d variants whose every expectation is quantitative, only %s of them non-disruptive, calling every "
          "variant destabilized on all three axes scores %s against %s for COMAVI, and the balanced score is %s against "
          "a ceiling of %s for any structure-blind strategy; of the %s grade points by which COMAVI exceeds calling no "
          "variant disrupted, the BRCT domain contributes %s and the hemoglobin tetramer %s (S1 Text). Localization also "
          "depends on how the %s expectations of a change measured only as an amount of complex are scored: without the "
          "evidence-resolution rule, %d rather than %d of the %d mechanisms are localized"
          % (len(qi), WORDS_[len(q_nd)], r3(q_all), r3(q_comavi), r3(q_bal), r3(ceiling, 1), r3(sum(contrib.values()), 1),
             r3(contrib["brca1_brct"], 1), r3(contrib["hemoglobin_tetramer"], 1), WORDS_[n_either], att_no_rule, att, n_m))
    # 5. resampling and permutation (independent seeds)
    g25 = np.array([float(GM[grade(label_axes(L), *p)]) for L, p in zip(B.mech_t25, PROFILES)])
    systems = B.system.values
    uniq = np.unique(systems)
    idx = {s: np.where(systems == s)[0] for s in uniq}
    rng = np.random.default_rng(20261112)
    boot = []
    for _ in range(args.draws):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        sel = np.concatenate([idx[s] for s in pick])
        boot.append(g25[sel].mean())
    lo, hi = np.percentile(boot, [2.5, 97.5])
    m = re.search(r"Resampling whole systems gave a 95% interval of (\d\.\d+)\u2013(\d\.\d+)", body)
    check("system-resampling interval (independent seed, %d draws)" % args.draws,
          m and abs(lo - float(m.group(1))) <= 0.012 and abs(hi - float(m.group(2))) <= 0.012,
          "%.3f-%.3f against the printed %s" % (lo, hi, m.group(0)[-11:] if m else "(not found)"))
    loso = [Fraction(int(round((g25.sum() - g25[idx[s]].sum()) * 2)), 2 * (len(B) - len(idx[s]))) for s in uniq]
    claim("leave one system out", "%s-%s" % (r3(min(loso)), r3(max(loso))),
          "leaving out one system at a time gave values from %s to %s" % (r3(min(loso)), r3(max(loso))))
    rng = np.random.default_rng(20261113)
    nulls = []
    labels25 = list(B.mech_t25)
    for _ in range(args.draws // 4):
        perm = rng.permutation(len(B))
        nulls.append(float(score(labels25, [PROFILES[i] for i in perm])) / len(B))
    nulls = np.array(nulls)
    m = re.search(r"permuting the ground truth gives a mean of (\d\.\d+) \(standard deviation (\d\.\d+); p = (\d\.\d+)", body)
    check("label-permutation null (independent seed, %d permutations)" % len(nulls),
          m and abs(nulls.mean() - float(m.group(1))) <= 0.004 and abs(nulls.std() - float(m.group(2))) <= 0.004
          and (nulls >= float(s25) / len(B) - 1e-12).sum() == 0,
          "mean %.3f, SD %.3f, %d permutations at or above the observed score; printed %s"
          % (nulls.mean(), nulls.std(), int((nulls >= float(s25) / len(B) - 1e-12).sum()), m.group(0)[-40:] if m else "(not found)"))
    # 6. direction agreement
    for t in TAGS:
        acc = collections.defaultdict(lambda: [0, 0])
        for (_, r), (exp, cls, either) in zip(B.iterrows(), PROFILES):
            for k, (a, n) in agreement(r, t, exp, either).items():
                acc[k][0] += a
                acc[k][1] += n
        k4, n4 = sum(v[0] for v in acc.values()), sum(v[1] for v in acc.values())
        stored = (int(B["structural_agreement_n_%s" % t].sum()), int(B["structural_agreement_d_%s" % t].sum()))
        check("four-output agreement at %s matches the stored per-variant counts" % t, (k4, n4) == stored,
              "%d/%d re-derived, %d/%d stored" % (k4, n4, *stored))
        if t == "t25":
            e_k = sum(acc[a][0] for a in ("monomer", "fold_complex", "binding", "either"))
            e_n = sum(acc[a][1] for a in ("monomer", "fold_complex", "binding", "either"))
            claim("direction agreement, isolated subunit", "%d/%d" % tuple(acc["monomer"]),
                  "%d of %d expectations for isolated-subunit stability" % tuple(acc["monomer"]))
            claim("direction agreement, assembled complex", "%d/%d" % tuple(acc["fold_complex"]),
                  "%d of %d for assembled-complex stability" % tuple(acc["fold_complex"]))
            claim("direction agreement, binding", "%d/%d" % tuple(acc["binding"]), "%d of %d for binding" % tuple(acc["binding"]))
            claim("direction agreement, joint row", "%d/%d" % tuple(acc["either"]), "met %d of them" % acc["either"][0])
            claim("direction agreement, energy axes", "%d/%d" % (e_k, e_n), "%d of the %d axes overall" % (e_k, e_n))
    # 7. priority score
    pv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
    pop = set(zip(pv.system, pv.variant))
    P = B[[k in pop for k in zip(B.system, B.variant)]].copy()
    EC = [priority(r) for _, r in P.iterrows()]
    P["E"], P["C"] = [e for e, _ in EC], [c for _, c in EC]
    P["S"] = (P.E + P.C) / 2
    stored_s = P.merge(pv[["system", "variant", "isds_v1"]].rename(columns={"isds_v1": "stored"}), on=["system", "variant"])
    check("priority score follows from its formula", np.allclose(stored_s.S, stored_s.stored, atol=1e-12),
          "max difference %.1e over %d variants" % (float((stored_s.S - stored_s.stored).abs().max()), len(P)))
    y = (P.expected_mech_class != SILENT).values
    order = np.argsort(-P.S.values, kind="stable")
    auc = roc_auc_score(y, P.S)
    claim("priority ROC AUC and average precision", r3(auc),
          "Its ROC AUC was %s and its average precision %s" % (r3(auc), r3(average_precision_score(y, P.S))))
    claim("top 10 and top 20", "%d, %d" % (y[order[:10]].sum(), y[order[:20]].sum()),
          "placed %d of the %d structural mechanisms in the top 10 and %d in the top 20"
          % (y[order[:10]].sum(), y.sum(), y[order[:20]].sum()))
    claim("energy component AUC", r3(roc_auc_score(y, P.E)), "The energy component alone reached an AUC of %s" % r3(roc_auc_score(y, P.E)))
    claim("context component AUC", r3(roc_auc_score(y, P.C)), "The context component alone reached %s" % r3(roc_auc_score(y, P.C)))
    # 8. tier screen
    strong = P.comavi_tier.isin(["Tier 1", "Tier 2"]).values
    tp, fp = int((strong & y).sum()), int((strong & ~y).sum())
    fn, tn = int(y.sum()) - tp, int((~y).sum()) - fp
    pf = fisher_exact([[tp, fn], [fp, tn]])[1]
    pp = exact_perm_p(P.system, y, strong)
    claim("tier screen", "%d/%d, %d/%d" % (tp, tp + fn, fp, fp + tn),
          "together with %d of the %d non-disruptive variants (specificity %s" % (fp, fp + tn, r3(Fraction(tn, fp + tn))))
    claim("tier screen p values", "%s, %s" % (r3(pf), r3(pp)), "Fisher p = %s, within-system permutation p = %s" % (r3(pf), r3(pp)))
    # 9. pathogenicity comparison
    cs = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_common_set.csv")[["system", "variant", "clinical_y"]]
    Q = P.merge(cs, on=["system", "variant"])  # P carries the canonical AlphaMissense column
    yq = Q.clinical_y.astype(int)
    claim("pathogenicity AUCs", ", ".join(r3(roc_auc_score(yq, Q[c])) for c in ("AM pathogenicity", "E", "S", "C")),
          "AlphaMissense ranked pathogenicity better (ROC AUC %s) than the energy component (%s), the priority score (%s) "
          "or the context component (%s)" % tuple(r3(roc_auc_score(yq, Q[c])) for c in ("AM pathogenicity", "E", "S", "C")))
    # 10. measured destabilizations recovered
    cal = pd.read_csv(RO / "COMAVI_delta_calibration_points.csv")
    d = cal[cal.measured_kcal >= 1.0]
    axis_i = {"Monomer fold": 0, "Multimer fold": 1, "Binding": 2}
    ax = d.axis.map(axis_i)
    assert ax.notna().all(), sorted(set(d.axis))
    rec = {t: int(sum(1 for v, i in zip(d.foldx_ddg, ax) if v >= THR[t][int(i)])) for t in TAGS}
    claim("measured destabilizations recovered", "%d, %d, %d of %d" % (rec["t10"], rec["t25"], rec["tSAP"], len(d)),
          "It did for %d at 1.0 kcal/mol, %d at 2.5 kcal/mol and %d at the prediction-interval limit"
          % (rec["t10"], rec["t25"], rec["tSAP"]))

    if args.report:
        pd.DataFrame(rows).to_csv(args.report, index=False)
    for f in fails:
        print("FAIL: " + f)
    if fails:
        return 1
    print("PASS: %d headline checks re-derived from the canonical table, the evidence ledger and the measured inputs, "
          "without importing any generator; every value is stated in the manuscript as re-derived" % len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
