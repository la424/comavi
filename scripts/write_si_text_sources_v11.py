#!/usr/bin/env python3
"""Write the v11 sources of S1-S3 Text from the committed records.

WHY THIS EXISTS
---------------
The v9.2 supporting texts were hand-edited Word files and kept values the
records had since replaced. This script writes their v11 sources the way the
v11 manuscript text was written: every value is read from a committed record,
every premise the prose relies on is asserted before the sentence is written,
and every value is listed with its record path in
docs/audit/COMAVI_v11_si_text_values.json. Rounding is done once, half to even,
on the exact value (grades are multiples of 0.5, so most statistics are exact
fractions); a value stored already rounded to one more place with a final 5 is
refused rather than rounded twice.

  records  ->  docs/manuscript/si/S{1,2,3}_Text_source.txt  ->  (build_si_texts.py)  ->  S{1,2,3}_Text.docx

Usage:
    python scripts/write_si_text_sources_v11.py           # write the sources and the values audit
    python scripts/write_si_text_sources_v11.py --check   # regenerate in memory, compare, write nothing
"""
import argparse
import collections
import importlib.util
import json
import pathlib
import re
import sys
from decimal import Decimal, ROUND_HALF_EVEN
from fractions import Fraction

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from sklearn.metrics import roc_auc_score

RD = pathlib.Path(__file__).resolve().parents[1]
RO = RD / "reference_outputs"
SI_SRC = RD / "docs" / "manuscript" / "si"
sys.path[:0] = [str(RD / "scripts"), str(RD)]
import apply_concordance_v5 as ac  # noqa: E402
import benchmark_population as bp  # noqa: E402
import build_detection_vs_attribution as bda  # noqa: E402
import build_manuscript as bm  # noqa: E402


def J(f):
    return json.loads((RO / f).read_text())


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


VALUES = {}


def keep(name, value, source):
    """Record a value with the record it came from; return it unchanged."""
    VALUES[name] = {"value": value if not isinstance(value, (np.integer, np.floating)) else value.item(),
                    "source": source}
    return value


# ------------------------------------------------------------------ rounding
MINUS = "\u2212"


class Ambiguous(ValueError):
    pass


def exact(x):
    """The exact rational behind a value. Grades are multiples of 0.5 and every
    reported statistic is a ratio or an interpolated order statistic of such
    means, so a float is snapped to the nearest fraction with a bounded
    denominator when it lies within float error of it. This matters only at
    rounding ties: the percentile bound -7/80 arrives as -0.08749999999999997,
    which rounds to -0.087, while the value it represents rounds half to even
    to -0.088."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, (int, np.integer)):
        return Fraction(int(x))
    f = Fraction(repr(float(x)))
    g = f.limit_denominator(10 ** 6)
    return g if abs(float(g) - float(x)) < 1e-12 else f


def r(x, nd, full=False):
    """Round once, half to even, on the exact value. A stored value already
    rounded to nd + 1 places whose last digit is 5 cannot be rounded correctly;
    pass the full-precision value (full=True) or a Fraction instead."""
    if not isinstance(x, Fraction) and not full:
        s = repr(float(x))
        if "e" not in s and "." in s:
            dec = s.split(".")[1]
            if len(dec) == nd + 1 and dec[-1] == "5":
                raise Ambiguous("%s rounded to %d places is ambiguous; supply full precision" % (s, nd))
    q = exact(x) if (full or isinstance(x, Fraction)) else Fraction(repr(float(x))) \
        if not isinstance(x, (int, np.integer)) else Fraction(int(x))
    neg = q < 0
    a = -q if neg else q
    scaled = a * 10 ** nd
    lo = scaled.numerator // scaled.denominator
    rem = scaled - lo
    if rem > Fraction(1, 2) or (rem == Fraction(1, 2) and lo % 2 == 1):
        lo += 1
    digits = str(lo).rjust(nd + 1, "0")
    out = digits[:-nd] + "." + digits[-nd:] if nd else digits
    if neg and lo != 0:
        out = MINUS + out
    return out


def ci(lo, hi, nd=3, full=False):
    """An interval; 'to' instead of an en dash when the lower bound is negative.
    full=True for bounds read at full precision, which are rounded half to even."""
    a, b = r(lo, nd, full), r(hi, nd, full)
    return "%s to %s" % (a, b) if a.startswith(MINUS) else "%s\u2013%s" % (a, b)


def frac(k, n, nd=3):
    kk = Fraction(str(k)) if isinstance(k, float) else Fraction(int(k))
    return "%s/%s = %s" % (num(k), num(n), r(kk / int(n), nd))


def num(k):
    if isinstance(k, float) and k.is_integer():
        return "%d" % k
    return ("%s" % k)


def half(x):
    """Grade totals are multiples of 0.5; print one decimal."""
    assert abs(x * 2 - round(x * 2)) < 1e-9, x
    return "%.1f" % x


def exact_mean(x, n):
    """A stored 4-decimal mean of n half-integer grades, as an exact Fraction."""
    tot = round(float(x) * n * 2) / 2
    assert abs(tot / n - float(x)) < 6e-5, (x, n)
    return Fraction(int(round(tot * 2)), 2 * n)


def exact_p(p, draws):
    """A stored fraction of resamples, as the exact count over the draws."""
    k = round(float(p) * draws)
    assert abs(k / draws - float(p)) < 5e-5, (p, draws)
    return Fraction(k, draws)


WORDS = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
         8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}
Cap = lambda s: s[0].upper() + s[1:]


def thousands(n):
    return "{:,}".format(int(n))


# ------------------------------------------------------------------ inputs
can = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
B = bp.benchmark(can).reset_index(drop=True)
assert len(B) == bp.N_BENCHMARK == 56
SIL = "structurally_silent"
mech = B.expected_mech_class != SIL
N, N_MECH, N_ND = len(B), int(mech.sum()), int((~mech).sum())
assert (N_MECH, N_ND, B.system.nunique()) == (29, 27, 14)
GM = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
TAGS = ["t10", "t15", "t20", "t25", "tSAP"]
TLAB = {"t10": "1.0", "t15": "1.5", "t20": "2.0", "t25": "2.5", "tSAP": "the prediction-interval limit"}
SYS = {sid: disp for disp, sid, _, _ in bm.SYSTEMS}
GENE = {r_.system + "|" + r_.variant: str(r_.gene).upper() for r_ in can.itertuples()}


def vn(system, variant):
    return "%s %s" % (GENE[system + "|" + variant], variant)


def grades(tag, mask=None):
    g = B["mech_consistency_%s" % tag].map(GM)
    assert g.notna().all()
    return g if mask is None else g[mask]


tot25 = grades("t25").sum()
assert tot25 == 37.0
keep("mechanism_pattern_total_t25", float(tot25), "scored_61var_canonical.csv mech_consistency_t25 over the 56 benchmark variants")
SCORE25 = Fraction(int(tot25 * 2), 2 * N)

refs = json.loads((RD / "docs/manuscript/references.json").read_text())
refs.update(json.loads((SI_SRC / "si_references.json").read_text()))


def check_keys(text):
    for m in bm.CITE.finditer(text):
        for k in bm.cite_keys(m.group(1)):
            assert k in refs, k


# ================================================================ S1 Text
# S1 Text: implementation details, statistical conventions and evidence provenance.
# Executed after v11_si_text_edits.py in the same namespace.

# ---- FoldX invocation: premises read from the shipped runner
runner = (RD / "scripts/comavi_v7/foldx_runner.py").read_text()
for needle in ('"--command=RepairPDB"', '"--command=BuildModel"', '"--command=AnalyseComplex"',
               'f"--analyseComplexChains={chains}"', '"--mutant-file=individual_list.txt"',
               'f"--numberOfRuns={n_runs}"', '_REPAIRED_CACHE'):
    assert needle in runner, needle
assert "AnalyzeComplex" not in runner          # v9.2 S1 Text misspelled the command
assert re.search(r"n_runs: int = 5", runner) or re.search(r"n_runs=5", runner), "five replicates"

# ---- tier constants: read from the shipped code, not typed
import comavi_v7.constants as cst  # noqa: E402
import comavi_v7.mechanism as mech_mod  # noqa: E402
cfg = load_module(RD / "scripts/comavi_v7/config.py", "comavi_cfg")
assert cfg.DISRUPTION_POINTS == [(20, 4.0), (10, 3.0), (4, 2.0), (1, 1.0)]
assert cfg.TIER_THRESHOLDS == [(5.0, "Tier 1"), (3.0, "Tier 2"), (1.5, "Tier 3")]
assert cst.grantham_severity(53.75) == 1.0 and cst.grantham_severity(1000) == 4.0
import comavi_v7.structure_loading as sload  # noqa: E402
assert [sload.classify_burial(x * cst.MAX_SASA["A"], "A") for x in (0.05, 0.10, 0.2, 0.25)] == \
    ["buried_core", "partially_buried", "partially_buried", "surface_exposed"]
msrc = (RD / "scripts/comavi_v7/mechanism.py").read_text()
for needle in ("if pl_plddt < 50:", "score += 2.0", "score += 1.5", "if n_gated >= 2:"):
    assert needle in msrc, needle
vt = (RD / "scripts/verify_tier_construction.py").read_text()
assert "score *= 0.4" in vt and "score *= 0.7" in vt and "< 50" in vt and "< 70" in vt
sl = (RD / "scripts/comavi_v7/structure_loading.py").read_text()
assert "def count_intra_contacts(structure, chain_id=\"A\", distance=5.0" in sl and "< 3:" in sl
assert "def count_interface(structure, my_chain, partner_chain, distance=5.0)" in sl
assert "ShrakeRupley" in sl

# ---- ledger composition on the benchmark
led = pd.read_csv(RO / "COMAVI_evidence_ledger.csv", keep_default_na=False, dtype=str)
bkeys = set(zip(B.system, B.variant))
LB = led[[k in bkeys for k in zip(led.system, led.variant)]]
LX = led[[k not in bkeys for k in zip(led.system, led.variant)]]
E = {"E1": "E1_quantitative_energetic", "E2": "E2_quantitative_functional", "E3": "E3_qualitative_experimental",
     "E4": "E4_population_frequency", "E5": "E5_inferred_no_axis_assay"}
typ = {k: int((LB.evidence_type == v).sum()) for k, v in E.items()}
ax = LB.axis.value_counts().to_dict()
assert len(LB) == 97 and len(led) == 99 and sum(typ.values()) == 97
assert set(LB[LB.evidence_type.isin([E["E4"], E["E5"]])].expected_token) == {"neutral"}
fc = LB[LB.axis == "fold_complex"]
fc_inf, fc_e5 = int((fc.evidence_directness == "inferred").sum()), int((fc.evidence_type == E["E5"]).sum())
n_var_axes = LB.groupby(["system", "variant"]).ngroups
excl_axes = sorted(vn(s, v) for s, v in zip(LX.system, LX.variant))
assert n_var_axes == N - 4 and excl_axes == ["BRCA1 R1699L", "BRCA1 R1699Q"]
for k, v in typ.items():
    keep("ledger_benchmark_%s" % k, v, "COMAVI_evidence_ledger.csv rows of the 56 benchmark variants")
keep("ledger_benchmark_axes", dict(ax), "COMAVI_evidence_ledger.csv")
keep("ledger_fold_complex_inferred", fc_inf, "COMAVI_evidence_ledger.csv evidence_directness")

# ---- direction agreement by evidence class (regenerated record, benchmark population)
st = J("COMAVI_evidence_stratified_agreement.json")
acc, crude, strat = st["accounting"], st["crude"], st["stratified"]["per_token"]
assert (acc["n_committed_axes"], acc["n_ddg_axes_gradeable"], acc["n_committed_but_gated_out"]) == (97, 79, 18)
assert crude["quant_ok"] + crude["soft_ok"] == 59 and crude["quant_n"] + crude["soft_n"] == 79
sa_mod = load_module(RD / "scripts/analyze_evidence_stratified_agreement.py", "strat_mod")
strata = [(strat[t]["quant_ok"], strat[t]["quant_n"] - strat[t]["quant_ok"],
           strat[t]["soft_ok"], strat[t]["soft_n"] - strat[t]["soft_ok"]) for t in ("destab", "neutral")]
or_mh, _, p_mh = sa_mod.mantel_haenszel(strata)
assert abs(or_mh - st["stratified"]["OR_mantel_haenszel"]) < 5e-4 and abs(p_mh - st["stratified"]["p"]) < 5e-4
_, p_crude = fisher_exact([[crude["quant_ok"], crude["quant_n"] - crude["quant_ok"]],
                           [crude["soft_ok"], crude["soft_n"] - crude["soft_ok"]]])
assert abs(p_crude - crude["p"]) < 5e-5
e1 = st["by_evidence_type"][E["E1"]]
keep("stratified_MH_OR_full_precision", or_mh, "recomputed from COMAVI_evidence_stratified_agreement.json strata "
     "with analyze_evidence_stratified_agreement.mantel_haenszel (stored value is rounded to 3 places)")
keep("crude_fisher_p_full_precision", p_crude, "fisher_exact on COMAVI_evidence_stratified_agreement.json crude counts")

# ---- tier screen on evidence-restricted subsets
tc = J("COMAVI_tier_construction.json")
v_all, v_e13, v_dc = (tc["evidence_restricted_views"][k] for k in ("all", "all_axes_E1_E3", "all_axes_direct_or_coupled"))
assert {k: v_e13[k] for k in ("n", "structural_n", "silent_n", "sensitivity", "specificity")} == \
       {k: v_dc[k] for k in ("n", "structural_n", "silent_n", "sensitivity", "specificity")}
rec_e13 = round(v_e13["sensitivity"] * v_e13["structural_n"])
assert Fraction(rec_e13, v_e13["structural_n"]) == Fraction(13, 16)
scr = tc["screen_full_tier"]
assert (scr["strong_structural"], scr["weak_structural"]) == (17, 4)

# ---- evidence-restricted views of the mechanism-pattern score
es = J("COMAVI_evidence_sensitivity.json")
assert (es["draws"], es["seed"], es["resampling_unit"], es["constant_labels_graded"]) == (2000, 8, "system", 31)
assert es["cross_check"]["stress_tests_system_cluster_ci95"] == es["cross_check"]["full_benchmark_score_ci95"]
VIEW_NAMES = [
    ("full_benchmark", "Full benchmark"),
    ("without_population_only_negatives", "Without non-disruptive variants supported only by population frequency"),
    ("without_whole_variant_negatives", "Without non-disruptive variants graded only at the whole-variant level"),
    ("without_weak_evidence_negatives", "Without all non-disruptive variants resting on population, inferred or whole-variant evidence"),
    ("every_expectation_E1_E3", "Every committed expectation E1\u2013E3"),
    ("every_expectation_E1_E2", "Every committed expectation E1\u2013E2"),
]
assert [k for k, _ in VIEW_NAMES] == list(es["views"])
NEUTRAL_CONST = {"No structural effect detected", "Interface variant (DDG neutral)",
                 "Structural variant \u2014 contact-driven (DDG neutral)", "Structural variant \u2014 burial-driven (DDG neutral)"}


def const_name(labels):
    labels = labels if isinstance(labels, list) else [labels]
    if set(labels) <= NEUTRAL_CONST:
        return "calling no variant disrupted (%d labels tie)" % len(labels) if len(labels) > 1 else "calling no variant disrupted"
    assert len(labels) == 1, labels
    return "\u201c%s\u201d for every variant" % labels[0]


# Intervals are re-rounded here to three decimals, so they are read at full
# precision from the generator itself (digits=None) and checked against the
# committed 4-decimal record, which they must reproduce.
bes = load_module(RD / "scripts/build_evidence_sensitivity.py", "bes_mod")
esF = bes.build(digits=None)


def _close(a, b):
    if isinstance(a, list):
        return len(a) == len(b) and all(_close(x, y) for x, y in zip(a, b))
    return abs(round(float(a), 4) - float(b)) < 1e-9


for key in es["views"]:
    for fld in ("comavi_score_ci95", "margin_ci95", "comavi_balanced", "comavi_balanced_ci95",
                "margin_p_at_or_below_zero"):
        assert _close(esF["views"][key][fld], es["views"][key][fld]), (key, fld)
view_rows = []
for key, label in VIEW_NAMES:
    v = dict(es["views"][key])
    for fld in ("comavi_score_ci95", "margin_ci95", "comavi_balanced", "comavi_balanced_ci95"):
        v[fld] = esF["views"][key][fld]
    assert v["best_constant_balanced"] <= 0.5 + 1e-9
    assert v["mechanisms"] + v["non_disruptive"] == v["n"]
    sc = exact_mean(v["comavi_score"], v["n"])
    stc = exact_mean(v["strongest_constant"], v["n"])
    mg = sc - stc
    assert abs(float(mg) - v["margin"]) < 6e-5
    view_rows.append([label, "%d (%d / %d)" % (v["n"], v["mechanisms"], v["non_disruptive"]), "%d" % v["systems"],
                      "%s (%s)" % (r(sc, 3), ci(*v["comavi_score_ci95"], full=True)),
                      "%s, %s" % (Cap(const_name(v["strongest_constant_labels"])), r(stc, 3)),
                      "%s (%s)" % (r(mg, 3), ci(*v["margin_ci95"], full=True)),
                      r(exact_p(v["margin_p_at_or_below_zero"], es["draws"]), 3),
                      "%s (%s)" % (r(v["comavi_balanced"], 3, full=True), ci(*v["comavi_balanced_ci95"], full=True))])
    keep("evidence_view_%s" % key, {k: v[k] for k in ("n", "comavi_score", "strongest_constant", "margin",
                                                      "margin_ci95", "margin_p_at_or_below_zero", "comavi_balanced")},
         "COMAVI_evidence_sensitivity.json views.%s" % key)
removed = {k: sorted(es["views"][k]["removed"]) if isinstance(es["views"][k]["removed"], list) else es["views"][k]["removed"]
           for k in ("without_population_only_negatives", "without_whole_variant_negatives")}

# ---- per-system contributions to the margin
bys = es["by_system"]
contrib = {s: d["contribution_vs_calling_nothing"] for s, d in bys.items()}
assert abs(sum(contrib.values()) - 10.0) < 1e-9
full = es["views"]["full_benchmark"]
assert abs(exact_mean(full["strongest_constant"], N) - Fraction(27, 56)) == 0
neg = sorted((s for s in contrib if contrib[s] < 0), key=lambda s: contrib[s])
zero = sorted(s for s in contrib if contrib[s] == 0)
false_calls = B[~mech & (B.mech_consistency_t25 == "inconsistent")]
fc_by_sys = false_calls.groupby("system").size().to_dict()
assert set(fc_by_sys) == set(neg), (fc_by_sys, neg)
assert all(-contrib[s] == fc_by_sys[s] - grades("t25", mech & (B.system == s)).sum() for s in neg)
loso_m = {s: esF["by_system"][s]["margin_without_system"] for s in bys}
assert all(_close(loso_m[s], bys[s]["margin_without_system"]) for s in bys)
lo_s, hi_s = min(loso_m, key=loso_m.get), max(loso_m, key=loso_m.get)
assert (lo_s, hi_s) == ("brca1_brct", "troponin_ic")
sys_rows = []
for s in sorted(bys, key=lambda s: (-contrib[s], SYS[s])):
    d = bys[s]
    sys_rows.append([SYS[s], "%d" % d["variants"], "%d" % d["non_disruptive"],
                     ("%+.1f" % d["contribution_vs_calling_nothing"]).replace("-", MINUS).replace("+0.0", "0.0"),
                     r(loso_m[s], 3, full=True)])
keep("per_system_contribution", contrib, "COMAVI_evidence_sensitivity.json by_system")

# ---- seeds and draw counts: read from the generators and their records
stress = load_module(RD / "verification/stress_tests.py", "stress_mod")
stsrc = (RD / "verification/stress_tests.py").read_text()
defaults = {m.group(1): int(m.group(2)) for m in re.finditer(r'add_argument\("--(n-perm|n-noise|n-boot)", type=int, default=(\d+)\)', stsrc)}
z = np.load(RO / "stress_tests/comavi_stress_draws.npz")
assert (len(z["mc_null"]), len(z["noise_mc"]), len(z["boot_mc"]), len(z["cluster_mc"])) == \
       (defaults["n-perm"], defaults["n-noise"], defaults["n-boot"], defaults["n-boot"])
assert int(z["seed_cluster"][0]) == stress.SEED_CLUSTER == es["seed"]
dva = J("COMAVI_detection_vs_attribution.json")
rs = J("COMAVI_rubric_sensitivity.json")
opd = J("COMAVI_operating_point_diagnostics.json")
sweep_mod = load_module(RD / "scripts/build_threshold_sweep.py", "sweep_mod")
cc = J("COMAVI_cohort_composition_checks.json")["structure_source"]["difference_experimental_minus_af3"]
sub_mod_src = (RD / "scripts/build_subthreshold_evidence.py").read_text()
sub_boot = int(re.search(r"^N_BOOT = (\d+)", sub_mod_src, re.M).group(1))
sub_seed = int(re.search(r"^SEED = (\d+)", sub_mod_src, re.M).group(1))
subrec = J("COMAVI_subthreshold_evidence.json")["axis_level"]
isb = J("isds_v1/ISDS_v1_summary.json")["bootstrap"]
seed_rows = [
    ["Mechanism-pattern score interval; margins, P(margin \u2264 0) and balanced scores of Table C", "Protein system",
     thousands(es["draws"]), "%d" % es["seed"], "build_evidence_sensitivity.py"],
    ["The same score interval, recomputed independently", "Protein system", thousands(len(z["cluster_mc"])),
     "%d" % stress.SEED_CLUSTER, "stress_tests.py"],
    ["Gain in detection rate, net gain and wrong-axis fraction", "Protein system",
     thousands(dva["system_cluster"]["draws"]), "%d" % dva["system_cluster"]["seed"], "build_detection_vs_attribution.py"],
    ["Half credit for calls on untested axes (S2 Text)", "Protein system", thousands(rs["draws"]), "%d" % rs["seed"],
     "build_rubric_sensitivity.py"],
    ["Specificity gain of the convergent evidence state (S2 Text)", "Protein system",
     thousands(tc["specificity_gain"]["cluster_bootstrap_n_draws"]), "0", "verify_tier_construction.py"],
    ["Paired operating-point comparison (S2 Text)", "Protein system", thousands(opd["paired_bootstrap"]["n_draws"]),
     "%d" % opd["paired_bootstrap"]["seed"], "build_operating_point_diagnostics.py"],
    ["Threshold-sweep bands (S3 Table)", "Protein system", thousands(sweep_mod.N_BOOT), "%d" % sweep_mod.SEED,
     "build_threshold_sweep.py"],
    ["Experimental against AlphaFold 3 systems", "Protein system, within source", thousands(cc["draws"]), "%d" % cc["seed"],
     "build_cohort_composition_checks.py"],
    ["Sub-threshold evidence AUC", "Protein system", "%s (%s with both classes present)" % (thousands(sub_boot), thousands(subrec["n_bootstrap_draws"])),
     "%d" % sub_seed, "build_subthreshold_evidence.py"],
    ["Priority-score AUC interval (S2 Text)", "Protein system", thousands(isb["draws_requested"]), "%d" % isb["seed"],
     "analyze_isds_v1.py"],
    ["Variant-level bootstrap intervals (S3 Text)", "Variant", thousands(len(z["boot_mc"])), "%d" % stress.SEED_BOOT, "stress_tests.py"],
    ["Label-permutation null (S2 Text)", "Expectation profile, permuted across variants", thousands(len(z["mc_null"])),
     "%d" % stress.SEED_PERM, "stress_tests.py"],
    ["Replicate-noise perturbation (S3 Text)", "Gaussian noise on each replicate-mean energy", thousands(len(z["noise_mc"])),
     "%d" % stress.SEED_NOISE, "stress_tests.py"],
    ["Stability of each call under noise (S3 Text)", "Gaussian noise on each replicate-mean energy", thousands(defaults["n-noise"]),
     "%d" % stress.SEED_FRAGILITY, "stress_tests.py"],
]
assert sub_seed == 11 and subrec["n_bootstrap_draws"] <= sub_boot

# ---------------------------------------------------------------- source text
def names_(xs):
    xs = list(xs)
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


def sys_list(keys):
    names = [SYS[k] for k in keys]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


S1 = []
A = S1.append
A("@title S1 Text. Implementation details, statistical conventions and evidence provenance.")
A("")
A("@contents")
A("")
A("## Section 1. FoldX calculations")
A("")
A("%from: scripts/comavi_v7/foldx_runner.py")
A("Each structure was repaired once with RepairPDB, and the repaired file was cached and reused for every variant "
  "scored on it. Variants were written in FoldX notation (wild-type residue, chain, residue number, substituted residue). "
  "BuildModel was run with five replicates on the isolated-subunit model and on the assembled complex. AnalyseComplex was "
  "run on the repaired wild-type complex and on each of the five mutant complexes, once for each pairing of the variant's "
  "chain with a partner chain. No FoldX option other than those shown was set, so every energy-function parameter is the "
  "FoldX 5.1 default [@schymkowitz2005; @delgado2019].")
A("")
A("@code foldx --command=RepairPDB --pdb=<structure> --output-dir=<workdir>")
A("@code foldx --command=BuildModel --pdb=<repaired structure> --mutant-file=individual_list.txt --numberOfRuns=5 --output-dir=<workdir>")
A("@code foldx --command=AnalyseComplex --pdb=<complex> --analyseComplexChains=<variant chain>,<partner chain> --output-dir=<workdir>")
A("")
A("%from: scripts/comavi_v7/foldx_runner.py")
A("Folding \u0394\u0394G is the mutant minus wild-type BuildModel energy averaged over the five replicates, so positive values "
  "are destabilizing. The wild-type interaction energy with a partner was calculated once, on the repaired complex, and the "
  "mutant interaction energy on each replicate mutant complex; binding \u0394\u0394G for that partner is the mean over "
  "replicates of mutant minus wild-type interaction energy. The replicate standard deviation is kept for every axis and partner.")
A("")
A("## Section 2. Structural-context tier: weights and cut-offs")
A("")
A("%from: scripts/comavi_v7/mechanism.py compute_disruption_score; scripts/comavi_v7/constants.py; scripts/comavi_v7/structure_loading.py; scripts/verify_tier_construction.py")
A("The tier is calculated from the structure alone and contains no energy. Substitution severity is the Grantham distance "
  "[@grantham1974] divided by 53.75 and capped at 4. The contact count of the variant residue is the number of residues of "
  "its own chain with any atom within 5 \u00c5 of any of its atoms, excluding the two sequence neighbours on either side, plus "
  "the largest number of residues within 5 \u00c5 in any one partner chain that passes the confidence gate. A partner passes the gate "
  "when the model confidence (pLDDT) at the variant site in that complex is at least 50; partners in experimental structures "
  "always pass. Burial uses the relative solvent accessibility of the residue, its Shrake\u2013Rupley accessible surface area "
  "[@shrake1973] calculated with Biopython [@cock2009] divided by the theoretical maximum for its residue type [@tien2013], "
  "and takes the most buried state across the isolated subunit and the complexes with confidence-passing partners. The best "
  "pLDDT at the site is the higher of its isolated-subunit and complex values. Table A gives the points and cut-offs.")
A("")
A("@tablecap Table A. Components and cut-offs of the structural-context tier.")
A("@tabhead Component | Definition | Points")
A("@tabrow Context | Severity \u00d7 contact count: 20 or more; 10 to <20; 4 to <10; 1 to <4; below 1 | 4; 3; 2; 1; 0")
A("@tabrow Interface | Contact with one confidence-passing partner; with two or more | 1.5; 2.0")
A("@tabrow Burial | Partially buried (relative accessibility 0.10 to <0.25); buried core (below 0.10) | 1.0; 2.0")
A("@tabrow Confidence | Best site pLDDT 50 to <70; below 50; experimental structure | \u00d70.7; \u00d70.4; \u00d71")
A("@tabrow Tier | Total of 5.0 or more; 3.0 to <5.0; 1.5 to <3.0; below 1.5 | Tier 1; Tier 2; Tier 3; Tier 4")
A("")
A("## Section 3. Statistical conventions, seeds and draw counts")
A("")
A("%from: scripts/build_evidence_sensitivity.py; verification/stress_tests.py; scripts/verify_tier_construction.py")
A("Proportions carry exact binomial 95% confidence intervals, and Fisher exact tests are two-sided. The within-system "
  "permutation tests of the tier screens are exact: each enumerates the convolution of the per-system hypergeometric "
  "distributions, so no seed is involved. Resampling intervals are percentile intervals. P(margin \u2264 0) is the fraction "
  "of resamples in which COMAVI's margin over the strongest structure-blind strategy is zero or negative, with that strategy "
  "re-selected among all 31 labels inside every resample. Values were computed at full precision and rounded once, half to "
  "even. Table B lists every resampling analysis with its unit, draw count and seed; scripts are in scripts/ except "
  "stress_tests.py, which is in verification/.")
A("")
A("@tablecap Table B. Resampling analyses, with unit, draw count and seed.")
A("@tabhead Analysis | Resampling unit | Draws | Seed | Script")
for row in seed_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 4. Evidence provenance and evidence-restricted analyses")
A("")
A("%from: reference_outputs/COMAVI_evidence_ledger.csv")
A("The 56 benchmark variants carry %d committed energetic-axis expectations (S1 Table): %d on isolated-subunit stability, "
  "%d on assembled-complex stability and %d on binding, on %d variants; the other four are graded only at the whole-variant "
  "level. By evidence class, %d are E1 quantitative energetic, %d E2 quantitative functional, %d E3 qualitative experimental, "
  "%d E4 population frequency and %d E5 inferred. The ledger also holds the two expectations of %s, which were scored but "
  "excluded from the benchmark (S2 Table, excluded_variants sheet). E4 and E5 commitments are heavily enriched for neutral "
  "expectations: all %d E4 axes and all %d E5 axes are neutral. Assembled-complex stability is the least independently "
  "grounded axis: %d of its %d committed axes rest on inferred rather than direct or coupled evidence, and %d are E5 inferred "
  "axis states."
  % (len(LB), ax["monomer"], ax["fold_complex"], ax["binding"], n_var_axes, typ["E1"], typ["E2"], typ["E3"], typ["E4"],
     typ["E5"], " and ".join(excl_axes), typ["E4"], typ["E5"], fc_inf, len(fc), fc_e5))
A("")
def _admitted(row, axis, need_interval):
    """Whether an expectation on this axis has a scoreable prediction: the
    monomer axis, or any partner on a complex axis, passes the confidence gate
    and, when need_interval, has a replicate interval excluding zero."""
    def ok(conf_col, dist_col, val_col):
        return (pd.notna(row.get(val_col)) and bool(row.get(conf_col)) is True
                and (not need_interval or bool(row.get(dist_col)) is True))
    if axis == "monomer":
        return ok("ddg_monomer_confident", "ddg_monomer_distinguishable_internal_from_0", "ddg_monomer")
    pre_ = "fold" if axis == "fold_complex" else "binding"
    return any(ok("ddg_%s_confident" % p, "ddg_%s_%s_distinguishable_internal_from_0" % (pre_, p), "ddg_%s_%s" % (pre_, p))
               for p in sorted({m.group(1) for c in B.columns
                                for m in [re.fullmatch(r"ddg_binding_([a-z0-9]+(?:_\d)?)", c)] if m}))


def _n_admitted(need_interval):
    n = 0
    for _, row in B.iterrows():
        exp = {r.axis: r.expected_token for r in LB[(LB.system == row.system) & (LB.variant == row.variant)].itertuples()}
        either = str(row.get("assembled_expectation_scope") or "") == ac.EITHER_ASSEMBLED_AXIS
        for a in exp:
            if either and a == "binding":
                n += _admitted(row, "binding", need_interval) or _admitted(row, "fold_complex", need_interval)
            else:
                n += _admitted(row, a, need_interval)
    return n


_n_conf, _n_both = _n_admitted(False), _n_admitted(True)
assert (_n_conf, _n_both) == (len(LB), acc["n_ddg_axes_gradeable"]), (_n_conf, _n_both)


A("%from: reference_outputs/COMAVI_evidence_stratified_agreement.json")
A("Of the %d expectations, %d can be compared with a prediction: on every energy axis they clear the model-confidence gate "
  "and the check that the interval of the predicted energy, the replicate mean \u00b1 1.96 replicate standard deviations, "
  "excludes zero. The other %d are removed by the interval check; the confidence gate removes none, and none is removed by "
  "missing curation. On these %d axes, direction agreement was "
  "%s for E1\u2013E2 evidence and %s for E3\u2013E5 evidence (Fisher p = %s). Evidence class was strongly associated with "
  "the expected axis state: direct measurements were preferentially available for destabilizers, whereas "
  "frequency-supported and inferred commitments were predominantly neutral. Within expected-destabilizing axes, agreement "
  "was %d/%d for E1\u2013E2 and %d/%d for E3\u2013E5; within expected-neutral axes, the corresponding rates were %d/%d and "
  "%d/%d. The Mantel\u2013Haenszel odds ratio across expected-state strata was %s (continuity-corrected test p = %s), so "
  "evidence class does not explain direction agreement once the expected axis state is considered. Among %d gradeable E1 "
  "energetic axes, %d agreed (%s)."
  % (acc["n_committed_axes"], acc["n_ddg_axes_gradeable"], acc["n_committed_but_gated_out"], acc["n_ddg_axes_gradeable"],
     frac(crude["quant_ok"], crude["quant_n"]), frac(crude["soft_ok"], crude["soft_n"]), r(p_crude, 3),
     strat["destab"]["quant_ok"], strat["destab"]["quant_n"], strat["destab"]["soft_ok"], strat["destab"]["soft_n"],
     strat["neutral"]["quant_ok"], strat["neutral"]["quant_n"], strat["neutral"]["soft_ok"], strat["neutral"]["soft_n"],
     r(or_mh, 2), r(p_mh, 2), e1["n"], e1["ok"], r(Fraction(e1["ok"], e1["n"]), 3)))
A("")
A("%from: reference_outputs/COMAVI_tier_construction.json evidence_restricted_views")
A("For the tier screen, requiring every committed axis to have E1\u2013E3 support leaves %d of the %d prioritization "
  "variants, %d structural and %d non-disruptive (tier sensitivity %s, 95%% CI %s; specificity %s, %s). Requiring every "
  "committed axis to be direct or coupled leaves the same %d variants and the same values. Recovery here rests on %d "
  "positives, of which %d are recovered, so the interval is wide; on the whole population the screen recovers %d of %d "
  "(sensitivity %s)."
  % (v_e13["n"], v_all["n"], v_e13["structural_n"], v_e13["silent_n"], r(Fraction(rec_e13, v_e13["structural_n"]), 3),
     ci(*v_e13["sensitivity_ci"]), r(v_e13["specificity"], 3), ci(*v_e13["specificity_ci"]), v_dc["n"],
     v_e13["structural_n"], rec_e13, scr["strong_structural"], scr["strong_structural"] + scr["weak_structural"],
     r(Fraction(scr["strong_structural"], scr["strong_structural"] + scr["weak_structural"]), 3)))
A("")
fb = es["views"]["full_benchmark"]
e12 = es["views"]["every_expectation_E1_E2"]
A("%from: reference_outputs/COMAVI_evidence_sensitivity.json views")
A("Table C reports the mechanism-pattern score on evidence-restricted subsets of the benchmark. Removing variants changes "
  "the balance of classes, and a structure-blind strategy scores, by construction, the fraction of a subset that it "
  "matches, so each subset is compared with its own strongest strategy, re-selected among the 31 labels on that subset. "
  "The balanced score, the mean of the mechanism and non-disruptive arms, does not depend on that balance: no structure-blind "
  "strategy exceeds 0.5 on it, because a label that accepts every non-disruptive variant recovers no mechanism, and a "
  "label that asserts disruption rejects no non-disruptive variant. On the %d variants whose expectations are all "
  "quantitative, only %s of them non-disruptive, calling every variant destabilized on all three axes outscores COMAVI, "
  "while COMAVI's balanced score stays above 0.5."
  % (e12["n"], WORDS[e12["non_disruptive"]]))
A("")
A("@tablecap Table C. Mechanism-pattern score on evidence-restricted subsets, each against its own strongest structure-blind strategy, at 2.5 kcal/mol.")
A("@tabhead Subset | Variants (mechanisms / non-disruptive) | Systems | COMAVI score (95% interval) | Strongest structure-blind strategy | Margin (95% interval) | P(margin \u2264 0) | Balanced score (95% interval)")
for row in view_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 5. Contribution of each system to the margin")
A("")
A("%from: reference_outputs/COMAVI_evidence_sensitivity.json by_system")
A("On the full benchmark the strongest structure-blind strategy is to call no variant disrupted, which scores %s of %d by "
  "accepting every non-disruptive variant. COMAVI scores %s, a margin of %s grade points. Table D divides that margin among "
  "the systems: a system's contribution is the sum, over its variants, of COMAVI's grade minus the grade of calling nothing. "
  "The BRCT domain contributes %s and the hemoglobin tetramer %s. Three systems reduce the margin, MSH2\u2013MSH6 and "
  "troponin I\u2013troponin C by %s each and PI3K by %s, because COMAVI predicts disruption for five non-disruptive variants "
  "there (%s) and earns no credit on the two PI3K mechanisms; %s systems contribute nothing. Leaving one system out, against "
  "the strongest strategy on the remaining variants, the margin ranges from %s without the BRCT domain to %s without "
  "troponin I\u2013troponin C."
  % (half(27.0), N, half(tot25), half(tot25 - 27.0), half(contrib["brca1_brct"]), half(contrib["hemoglobin_tetramer"]),
     half(-contrib["msh2_msh6"]), half(-contrib["pi3k"]), names_(sorted(vn(s, v) for s, v in zip(false_calls.system, false_calls.variant))),
     WORDS[len(zero)], r(loso_m[lo_s], 3, full=True), r(loso_m[hi_s], 3, full=True)))
assert sorted(neg) == ["msh2_msh6", "pi3k", "troponin_ic"]
assert contrib["msh2_msh6"] == contrib["troponin_ic"] == -2.0 and contrib["pi3k"] == -1.0
assert grades("t25", mech & (B.system == "pi3k")).sum() == 0 and int((mech & (B.system == "pi3k")).sum()) == 2
A("")
A("@tablecap Table D. Contribution of each system to COMAVI's margin over calling no variant disrupted, at 2.5 kcal/mol.")
A("@tabhead System | Variants | Non-disruptive | Contribution (grade points) | Margin without the system")
for row in sys_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("@references")
SOURCES = {"S1_Text": "\n".join(S1) + "\n"}


# ================================================================ S2 Text
# S2 Text: sensitivity analyses, mechanism-class decomposition and operating-point diagnostics.
# Executed after v11_si_text_edits.py and v11_si_s1.py in the same namespace.

CLASS_NAMES = [("fold_mechanism", "Fold"), ("mixed_structural", "Mixed fold and interface"),
               ("ppi_destab_mechanism", "Interface destabilization"), ("ppi_stab_mechanism", "Interface stabilization"),
               ("structurally_silent", "Non-disruptive")]

# ---- Section 1: structure-blind strategies and the permutation null
rb = J("COMAVI_rubric_baselines.json")
ecl = rb["every_constant_label"]
labs = ecl["labels"]
assert ecl["n_labels_graded"] == len(labs) == 31
# Each label's score must be its grade counts, exactly.
for l, v in labs.items():
    gc = v["grade_counts"]
    assert sum(gc.values()) == N
    assert abs((gc.get("consistent", 0) + 0.5 * gc.get("partial", 0)) / N - v["score"]) < 6e-5, l
named = {"never_fire": "No structural effect detected", "always_fold_and_ppi": "Both fold + PPI destabilization"}
for k, l in named.items():
    assert abs(rb["baselines"][k]["score"] - labs[l]["score"]) < 1e-9, k
lab_total = {l: v["grade_counts"].get("consistent", 0) + 0.5 * v["grade_counts"].get("partial", 0) for l, v in labs.items()}
best_total = max(lab_total.values())
ties = sorted(l for l, t in lab_total.items() if t == best_total)
assert best_total == 27.0 and set(ties) == NEUTRAL_CONST
assert lab_total["Both fold + PPI destabilization"] == 24.5
assert ecl["best_balanced"] <= 0.5 + 1e-9
fbv = es["views"]["full_benchmark"]
arm_m, arm_s = Fraction(int(round(fbv["arm_mechanism"] * N_MECH * 2)), 2 * N_MECH), Fraction(int(round(fbv["arm_non_disruptive"] * N_ND * 2)), 2 * N_ND)
assert abs(float(arm_m) - fbv["arm_mechanism"]) < 6e-5 and abs(float(arm_s) - fbv["arm_non_disruptive"]) < 6e-5
assert grades("t25", mech).sum() == float(arm_m * N_MECH) and grades("t25", ~mech).sum() == float(arm_s * N_ND)

const_rows = []
for l in sorted(labs, key=lambda l: (-lab_total[l], l)):
    v = labs[l]
    gc = v["grade_counts"]
    am = Fraction(int(round(v["arm_committed_structural"] * N_MECH * 2)), 2 * N_MECH)
    asl = Fraction(int(round(v["arm_curated_silent"] * N_ND * 2)), 2 * N_ND)
    assert abs(float(am) - v["arm_committed_structural"]) < 6e-5 and abs(float(asl) - v["arm_curated_silent"]) < 6e-5
    const_rows.append([l, r(Fraction(int(lab_total[l] * 2), 2 * N), 3), r(am, 3), r(asl, 3), r((am + asl) / 2, 3),
                       "%d / %d / %d" % (gc.get("consistent", 0), gc.get("partial", 0), gc.get("inconsistent", 0))])
keep("constant_labels", {l: lab_total[l] for l in labs}, "COMAVI_rubric_baselines.json every_constant_label.labels grade_counts")

n_assert = int((~B.mech_t25.isin(list(NEUTRAL_CONST))).sum())
assert n_assert == det_k_plus_false if "det_k_plus_false" in dir() else True
mc_null = z["mc_null"]
mc_obs = float(z["mc_obs"][0])
assert abs(mc_obs - float(SCORE25)) < 1e-12 and len(mc_null) == 2000
n_ge = int((mc_null >= mc_obs - 1e-12).sum())
p_perm = Fraction(n_ge + 1, len(mc_null) + 1)
st_csv = pd.read_csv(RO / "stress_tests/comavi_stress_tests.csv").set_index("test")
assert abs(float(st_csv.loc["permutation null (MC)", "p_value"]) - float(p_perm)) < 5e-5
null_mean, null_sd = float(mc_null.mean()), float(mc_null.std())
assert abs(round(null_mean, 4) - float(st_csv.loc["permutation null (MC)", "null_mean"])) < 1e-9
assert abs(round(null_sd, 4) - float(st_csv.loc["permutation null (MC)", "null_sd"])) < 1e-9
keep("permutation_null", {"mean": null_mean, "sd": null_sd, "n_ge_observed": n_ge, "p": float(p_perm)},
     "reference_outputs/stress_tests/comavi_stress_draws.npz mc_null (full precision)")

# ---- Section 2: grades by mechanism class
pc = dva["per_class"]
assert sum(pc[c]["n"] for c, _ in CLASS_NAMES) == N
for c, _ in CLASS_NAMES:
    sub = B[B.expected_mech_class == c]
    cnt = sub.mech_consistency_t25.value_counts().to_dict()
    assert (cnt.get("consistent", 0), cnt.get("partial", 0), cnt.get("inconsistent", 0)) == \
           (pc[c]["consistent"], pc[c]["partial"], pc[c]["inconsistent"]), c
class_rows = []
for c, name in CLASS_NAMES:
    v = pc[c]
    tot = v["consistent"] + 0.5 * v["partial"]
    det = ("%d" % v["detection"]) if c != SIL else "\u2013"
    loc = ("%d" % v["attribution"]) if c != SIL else "\u2013"
    rej = "\u2013" if c != SIL else "%d" % v["correct_rejection"]
    class_rows.append([name, "%d" % v["n"], "%d / %d / %d" % (v["consistent"], v["partial"], v["inconsistent"]),
                       "%s (%s)" % (half(tot), r(Fraction(int(tot * 2), 2 * v["n"]), 3)), det, loc, rej])
icd = dva["interface_class_diagnosis"]
assert (icd["n"], icd["binding_axis_fires"], icd["complex_fold_axis_fires"]) == (11, 4, 7)
fold_bad = B[(B.expected_mech_class == "fold_mechanism") & (B.mech_consistency_t25 != "consistent")]
fold_missed = fold_bad[fold_bad.mech_t25 == "No structural effect detected"]
assert sorted(vn(s, v) for s, v in zip(fold_missed.system, fold_missed.variant)) == \
       ["BRCA1 L1664P", "BRCA1 R1751Q", "BRCA1 V1736A", "MLH1 H718Y"]

# ---- Section 3: localization, axis designs and the evidence-resolution rule
det, att = dva["detection"], dva["attribution"]
assert (det["k"], det["n"], att["k"]) == (20, 29, 18)
msd = dva["missed_at_reference"]
assert msd["n_missed"] == det["n"] - det["k"] == 9
missed_names = ["%s %s" % (m["gene"], m["variant"]) for m in msd["variants"]]
axd = dva["axis_designs"]
design_rows = []
for key, name in (("isolated_subunit", "Isolated subunit only"),
                  ("isolated_subunit_and_binding", "Isolated subunit and binding"),
                  ("isolated_subunit_and_assembled_stability", "Isolated subunit and assembled-complex stability"),
                  ("all_three", "All three axes (COMAVI)")):
    d = axd[key]
    design_rows.append([name, "%d of %d" % (d["detection"], d["n_structural"]), "%d of %d" % (d["attribution"], d["n_structural"]),
                        "%d of %d" % (d["false_calls"], d["n_silent"])])
only_asm = sorted("%s %s" % (d["gene"], d["variant"]) for d in axd["detected_only_through_assembled_stability"])
sc_ = dva["system_cluster"]
# evidence-resolution rule
er = pd.read_csv(RD / "docs/audit/COMAVI_v11_evidence_resolution.csv")
erp = er[er.pipeline == "primary"]
assert set(erp.old_grade) == {"inconsistent"} and set(erp.new_grade) == {"consistent"}
gain_by_t = {t: int((erp.threshold == t).sum()) for t in TAGS}
assert gain_by_t == {"t10": 2, "t15": 1, "t20": 2, "t25": 3, "tSAP": 3}
assert sorted(erp[erp.threshold == "t25"].variant) == ["D351H", "L749P", "V11G"]
scope = B[B.assembled_expectation_scope.fillna("") == ac.EITHER_ASSEMBLED_AXIS]
scope_names = sorted(vn(s, v) for s, v in zip(scope.system, scope.variant))
assert len(scope_names) == 6
assay = {vn(s, v): a for s, v, a in zip(LB.system, LB.variant, LB.complex_amount_assay) if a}
neutral_amount = sorted(k for k in assay if k not in scope_names)
assert len(assay) == 10 and len(neutral_amount) == 4
_tok = set(t.strip() for a in assay.values() for t in a.split(";"))
assert _tok == {"Y2H", "mammalian two-hybrid", "co-IP", "MAVE interaction score", "MSH6 IP", "GST pull-down", "TR-FRET"}, _tok

fr_ = bda.build(can.copy(), return_frame=True)
STRUCT_CLASSES = ("fold_mechanism", "mixed_structural", "ppi_destab_mechanism", "ppi_stab_mechanism")
_F = fr_[1][fr_[1].expected_mech_class.isin(STRUCT_CLASSES)]
_dn = sorted(vn(s, v) for s, v, a, rt in zip(_F.system, _F.variant, _F.any_fires, _F.right_axis) if a and rt is not True and not (type(rt).__name__ == "bool_" and bool(rt)))
assert _dn == ["MLH1 R755W", "PIK3CA E542K"], _dn
assert abs(dva["attribution"]["rate"] - round(att["k"] / att["n"], 3)) < 1e-9
can_norule = can.copy()
can_norule["assembled_expectation_scope"] = ""
att0 = bda.build(can_norule)["attribution"]
assert (att0["k"], att0["n"]) == (15, 29)
sw = {row["threshold_tag"]: row for row in J("COMAVI_threshold_sweep.json")["rows"]}
tot_by_t = {t: grades(t).sum() for t in TAGS}
for t in TAGS:
    assert abs(sw[t]["whole_variant"] - tot_by_t[t] / N) < 6e-5
asm_either = (sw["t25"]["axis_assembled_either_k"], sw["t25"]["axis_assembled_either_n"])
assert asm_either == (4, 6)
keep("evidence_resolution_effect", {"grade_points_gained": gain_by_t, "attribution_without_rule": att0["k"]},
     "docs/audit/COMAVI_v11_evidence_resolution.csv; build_detection_vs_attribution.build without the scope column")

# ---- Section 4: half credit for calls on untested axes
rsr = rs["reference_2_5"]
gch = pd.DataFrame(rs["grade_changes"])
assert set(gch.primary_grade) == {"inconsistent"} and set(gch.sensitivity_grade) == {"partial"}
assert abs(rsr["score"] - (tot25 + 0.5) / N) < 6e-5 and abs(rsr["margin"] - (tot25 + 0.5 - 27) / N) < 1e-9
by_t_rs = rs["by_threshold"]
lab_t = {"1.0": "t10", "1.5": "t15", "2.0": "t20", "2.5": "t25", "2.9/2.9/3.5": "tSAP"}
for k, v in by_t_rs.items():
    assert abs(v["primary_score"] - tot_by_t[lab_t[k]] / N) < 6e-5, k
chg_by = gch.groupby("threshold").apply(lambda d: sorted("%s %s" % (g, v) for g, v in zip(d.gene, d.variant))).to_dict()

# ---- Section 5: operating points
cr = opd["class_reweighting"]
pb = opd["paired_bootstrap"]
assert (pb["n_draws"], pb["seed"]) == (10000, 11) and opd["population"]["n_graded"] == N
for t in TAGS:
    pe = opd["point_estimates"][t]
    assert abs(pe["recovery"] - grades(t, mech).mean()) < 6e-5, t
    assert abs(pe["correct_rejection"] - grades(t, ~mech).mean()) < 6e-5, t


def pct_of(p):
    q = Fraction(str(p))
    return "every resample" if q == 1 else "%s%% of resamples" % r(q * 100, 1)


op_rows = []
for t in TAGS:
    pe = opd["point_estimates"][t]
    op_rows.append([("%s kcal/mol" % TLAB[t]) if t != "tSAP" else "Prediction-interval limit (2.9 / 2.9 / 3.5)",
                    "%s (%s)" % (r(exact_mean(pe["recovery"], N_MECH), 3), ci(*pb["recovery"][t]["ci95"])),
                    "%s (%s)" % (r(exact_mean(pe["correct_rejection"], N_ND), 3), ci(*pb["correct_rejection"][t]["ci95"])),
                    "%s (%s)" % (half(tot_by_t[t]), r(Fraction(int(tot_by_t[t] * 2), 2 * N), 3))])

# ---- Section 6: context comparators and interface ablation (46-variant population)
cbt = pd.read_csv(RO / "isds_v1/ISDS_v1_component_binary_baselines.csv").set_index("screen")
assert set(cbt.within_system_permutation_method) == {"exact (convolved hypergeometrics)"}
full_s, nob = tc["screen_full_tier"], tc["screen_interface_bonus_zeroed"]
assert (cbt.loc["full_tier_1_2", "tp"], cbt.loc["full_tier_1_2", "tn"]) == (full_s["strong_structural"], full_s["weak_silent"])
assert (cbt.loc["tier_without_interface_bonus", "tp"], cbt.loc["tier_without_interface_bonus", "tn"]) == \
       (nob["strong_structural"], nob["weak_silent"])
_, fisher_iface = fisher_exact([[cbt.loc["interface_status_alone", "tp"], cbt.loc["interface_status_alone", "fn"]],
                                [cbt.loc["interface_status_alone", "fp"], cbt.loc["interface_status_alone", "tn"]]])
from math import comb as _comb


def exact_perm_p(strong):
    """Exact within-system permutation p as a Fraction (the method of
    scripts/verify_tier_construction.exact_within_system_perm)."""
    pop = pd.DataFrame({"system": pvv.system, "truth": pvv.structural_ground_truth.astype(bool), "strong": strong.astype(bool)})
    obs = int((pop["strong"] & pop["truth"]).sum())
    tot = {0: Fraction(1)}
    for _, g in pop.groupby("system"):
        n_, k_, ns_ = len(g), int(g["truth"].sum()), int(g["strong"].sum())
        d = {x: Fraction(_comb(ns_, x) * _comb(n_ - ns_, k_ - x), _comb(n_, k_))
             for x in range(max(0, k_ - (n_ - ns_)), min(k_, ns_) + 1)}
        nt = {}
        for s0, p0 in tot.items():
            for x, px in d.items():
                nt[s0 + x] = nt.get(s0 + x, Fraction(0)) + p0 * px
        tot = nt
    return sum(p_ for s_, p_ in tot.items() if s_ >= obs), obs


pvv = pd.read_csv(RO / "isds_v1/ISDS_v1_per_variant.csv")
screens = {"interface_status_alone": pvv.interface_only_positive.astype(bool),
           "tier_without_interface_bonus": pvv.no_interface_score >= 3.0,
           "full_tier_1_2": pvv.reported_tier.isin(["Tier 1", "Tier 2"])}
perm_p = {}
for key, strong in screens.items():
    pp_, obs_ = exact_perm_p(strong)
    q = cbt.loc[key]
    assert obs_ == q.tp and int(strong.sum()) == q.tp + q.fp, key
    assert abs(round(float(pp_), 4) - q.within_system_permutation_p) < 1e-9, (key, float(pp_))
    perm_p[key] = pp_
keep("within_system_permutation_p", {k: float(v) for k, v in perm_p.items()}, "exact convolution, recomputed here")
cmp_rows = []
def _fisher(key):
    q = cbt.loc[key]
    return fisher_exact([[int(q.tp), int(q.fn)], [int(q.fp), int(q.tn)]])[1]


assert abs(_fisher("full_tier_1_2") - full_s["fisher_p"]) < 5e-6 and abs(_fisher("tier_without_interface_bonus") - nob["fisher_p"]) < 5e-6
for key, name in (("interface_status_alone", "Interface status alone"),
                  ("tier_without_interface_bonus", "Tier 1\u20132 without the interface bonus"),
                  ("full_tier_1_2", "Tier 1\u20132 (full tier)")):
    fisher_p = _fisher(key)
    q = cbt.loc[key]
    cmp_rows.append([name, "%d of %d (%s)" % (q.tp, q.tp + q.fn, r(Fraction(int(q.tp), int(q.tp + q.fn)), 3)),
                     "%d of %d (%s)" % (q.tn, q.tn + q.fp, r(Fraction(int(q.tn), int(q.tn + q.fp)), 3)),
                     r(fisher_p, 3, full=True), r(perm_p[key], 3)])
reach = tc["interface_bonus_reach"]
assert (reach["structural_with_bonus"], reach["structural_total"]) == (15, 21)
assert len(pvv) == 46 and pvv.system.nunique() == 13
tier_frac = pvv.groupby("reported_tier").structural_ground_truth.agg(["sum", "count"])
tf = {t: (int(tier_frac.loc[t, "sum"]), int(tier_frac.loc[t, "count"])) for t in ("Tier 1", "Tier 2", "Tier 3", "Tier 4")}
assert tf == {"Tier 1": (11, 13), "Tier 2": (6, 16), "Tier 3": (4, 10), "Tier 4": (0, 7)}
dem = sorted("%s %s" % (GENE[d["system"] + "|" + d["variant"]], d["variant"]) for d in tc["all_demotions"])
assert all(d["comavi_tier"] == "Tier 2" and d["tier_no_iface"] == "Tier 3" for d in tc["all_demotions"])
no_bonus = sorted(vn(s, v) for s, v in zip(pvv.system, pvv.variant) if v in reach["structural_without_bonus"])

# ---- Section 7: evidence states
fes = J("isds_v1/ISDS_v1_evidence_state_tables.json")
states = ["Convergent", "Energy only", "Context only", "Neither"]
state_rows = []
for row in fes["four_evidence_states"]:
    t = row["threshold"]
    assert sum(row[s]["structural"] for s in states) == 21 and sum(row[s]["no_lesion"] for s in states) == 25
    state_rows.append([("%s kcal/mol" % TLAB[t]) if t != "tSAP" else "Prediction-interval limit"] +
                      ["%d / %d" % (row[s]["structural"], row[s]["no_lesion"]) for s in states])
dm25 = [d for d in fes["detection_metrics"] if d["threshold"] == "t25"][0]
a_, c_ = dm25["any_energetic_firing"], dm25["convergent_energy_and_tier"]
assert a_["tp"] == c_["tp"] == 15 and (a_["tn"], c_["tn"]) == (20, 23)
sg = tc["specificity_gain"]
assert sg["cluster_bootstrap_n_draws"] == 200000 and sg["by_threshold"]["t25"]["energy_only_variants"] == ["A636P", "R145G", "R145Q"]
eo = sorted(vn(s, v) for s, v in zip(B.system, B.variant) if v in sg["by_threshold"]["t25"]["energy_only_variants"])

# ---- Section 8: priority-score formulas, intervals and bands
psum = J("COMAVI_priority_sensitivity_summary.json")
y46 = pvv.structural_ground_truth.astype(bool)
npairs = int(y46.sum() * (~y46).sum())
assert npairs == 525
def pairs(col):
    a = roc_auc_score(y46, pvv[col]) * npairs
    assert abs(a * 2 - round(a * 2)) < 1e-9
    return round(a * 2) / 2
rep_pairs = pairs("isds_v1")
assert rep_pairs == 438.0 and abs(rep_pairs / npairs - psum["reported_roc_auc"]) < 6e-5
alt_rows, above, equal, rho_ctx = [], [], [], []
for row in psum["alternative_formulas"]["rows"]:
    pp = pairs(row["column"])
    assert abs(pp / npairs - row["roc_auc"]) < 6e-5
    if pp > rep_pairs: above.append((row["alternative"], pp))
    if pp == rep_pairs: equal.append(row["alternative"])
    rho = spearmanr(pvv["isds_v1"], pvv[row["column"]])[0]
    assert abs(rho - row["spearman_vs_reported"]) < 6e-5
    if row["graded_context_term"]:
        rho_ctx.append(rho)
    alt_rows.append([Cap(row["alternative"]), "yes" if row["graded_context_term"] else "no",
                     "%s (%s of %d pairs)" % (r(Fraction(int(pp * 2), 2 * npairs), 3), num(pp), npairs), r(rho, 3)])
assert len(above) == 5 and equal == ["cohort-standardized"]
top_alt = max(above, key=lambda x: x[1])
isb = J("isds_v1/ISDS_v1_summary.json")["bootstrap"]
loso = pd.DataFrame(psum["leave_one_system_out"]["folds"])
assert loso.combined_above_both.all() and len(loso) == 13
loso_exact = {}
for s_ in loso.dropped_system:
    keep_ = pvv.system != s_
    yy = y46[keep_]
    npp = int(yy.sum() * (~yy).sum())
    a2 = roc_auc_score(yy, pvv.loc[keep_, "isds_v1"]) * npp * 2
    assert abs(a2 - round(a2)) < 1e-9
    loso_exact[s_] = Fraction(int(round(a2)), 2 * npp)
    assert abs(float(loso_exact[s_]) - float(loso.set_index("dropped_system").combined[s_])) < 6e-5, s_
bands = pd.read_csv(RO / "isds_v1/ISDS_v1_score_bands.csv")
assert bands.n_variants.sum() == 46 and bands.n_structural_mechanism.sum() == 21
band_rows = [[b.isds_band.replace("-", "\u2013"), "%d" % b.n_variants, "%d" % b.n_structural_mechanism,
              r(Fraction(int(b.n_structural_mechanism), int(b.n_variants)), 3)] for b in bands.itertuples()]

# ---------------------------------------------------------------- source text
def names(xs):
    xs = list(xs)
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


S2 = []
A = S2.append
A("@title S2 Text. Sensitivity analyses, mechanism-class decomposition and operating-point diagnostics.")
A("")
A("@contents")
A("")
A("## Section 1. Structure-blind strategies and the permutation null")
A("")
A("%from: reference_outputs/COMAVI_rubric_baselines.json; reference_outputs/stress_tests")
A("A structure-blind strategy assigns the same label to every variant, without looking at any structure, and is graded "
  "with the same rubric as COMAVI (Table 3). Table A grades each of the 31 labels COMAVI can assign. Calling no variant "
  "disrupted scores %s of %d (%s): it accepts all %d non-disruptive variants and recovers none of the %d mechanisms. The three "
  "labels that name a structural or interface site without an energy change are graded as no-effect calls and score the "
  "same. The strongest label that asserts disruption, a fold and binding destabilization on every variant, scores %s (%s). "
  "COMAVI scores %s (%s), %s on the mechanisms and %s on the non-disruptive variants. Because a constant label either "
  "accepts every non-disruptive variant and recovers no mechanism or does the reverse, no constant label exceeds 0.5 on "
  "the balanced score, the mean of the two arms; COMAVI's balanced score is %s."
  % (half(best_total), N, r(Fraction(int(best_total * 2), 2 * N), 3), N_ND, N_MECH, half(lab_total["Both fold + PPI destabilization"]),
     r(Fraction(int(lab_total["Both fold + PPI destabilization"] * 2), 2 * N), 3), half(tot25), r(SCORE25, 3), r(arm_m, 3), r(arm_s, 3),
     r((arm_m + arm_s) / 2, 3)))
A("")
A("%from: verification/stress_tests.py; reference_outputs/stress_tests/comavi_stress_draws.npz")
A("The permutation null keeps every COMAVI call fixed and reassigns whole expectation profiles, the curated class with its "
  "per-axis expectations and their evidence annotations, among the %d variants. Over %s permutations the score had a mean "
  "of %s (standard deviation %s), and %s reached the observed %s, so p = %s, the smallest value %s permutations can give. "
  "The null mean lies below the strongest structure-blind strategy because a permutation keeps COMAVI's calls, which assert "
  "disruption for %d variants, while a constant label can exploit the class balance of the benchmark. The margin over the "
  "strongest strategy (S1 Text) is therefore the stricter comparison."
  % (N, thousands(len(mc_null)), r(null_mean, 3, full=True), r(null_sd, 3, full=True), "none" if n_ge == 0 else "%d" % n_ge,
     r(SCORE25, 3), r(p_perm, 4), thousands(len(mc_null)), n_assert))
A("")
A("@tablecap Table A. Every structure-blind strategy on the 56 benchmark variants at 2.5 kcal/mol, ordered by score.")
A("@tabhead Label assigned to every variant | Score | Mechanism arm (29) | Non-disruptive arm (27) | Balanced score | Consistent / partial / inconsistent")
for row in const_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 2. Grades by mechanism class")
A("")
A("%from: reference_outputs/COMAVI_detection_vs_attribution.json per_class")
A("Table B divides the grades by curated mechanism class. The non-disruptive class scores highest and the interface classes "
  "lowest. Among the %d variants curated as interface destabilizers, the binding axis reached the threshold for %d and the "
  "assembled-complex stability axis for %d, so these disruptions are more often detected than placed on the binding axis. "
  "Of the %d fold mechanisms not graded consistent, %d produce no call at all (%s); their measured or reported effects are "
  "described in the main text and S3 Table (missed_mechanisms sheet)."
  % (icd["n"], icd["binding_axis_fires"], icd["complex_fold_axis_fires"], len(fold_bad), len(fold_missed),
     names(sorted(vn(s, v) for s, v in zip(fold_missed.system, fold_missed.variant)))))
A("")
A("@tablecap Table B. Grades, detection and localization by curated mechanism class at 2.5 kcal/mol.")
A("@tabhead Class | Variants | Consistent / partial / inconsistent | Grade total (mean) | Detected | Localized | Correctly rejected")
for row in class_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 3. Localization, axis design and the evidence-resolution rule")
A("")
A("%from: reference_outputs/COMAVI_detection_vs_attribution.json axis_designs, system_cluster")
A("A mechanism is detected when any axis reaches the threshold, and localized when an axis the literature implicates "
  "reaches it; whether the predicted change has the expected direction is assessed by the grade. At 2.5 kcal/mol COMAVI "
  "detected %d of the %d mechanisms and localized %d. The two detected but not localized are PIK3CA E542K, an interface "
  "destabilizer that reaches the threshold only on assembled-complex stability, and MLH1 R755W, a fold mechanism that "
  "reaches it only on binding, where it is called as strengthened binding. Table C repeats the analysis with axes removed. The isolated subunit alone detects %d mechanisms; "
  "the three axes together detect %d, at a cost of %d further false calls among the %d non-disruptive variants. Five "
  "mechanisms are detected only through assembled-complex stability: %s. Resampling whole systems puts a 95%% interval of "
  "%s on the gain in detection rate, %s on the gain net of the added false calls and %s on the wrong-axis fraction."
  % (det["k"], det["n"], att["k"], axd["isolated_subunit"]["detection"], axd["all_three"]["detection"],
     axd["all_three"]["false_calls"] - axd["isolated_subunit"]["false_calls"], N_ND, names(only_asm),
     ci(*sc_["detection_gain"]["ci95"], nd=2), ci(*sc_["net_gain"]["ci95"], nd=2), ci(*sc_["localization_gap"]["ci95"], nd=2)))
A("")
A("@tablecap Table C. Detection, localization and false calls by axis design at 2.5 kcal/mol.")
A("@tabhead Axes used | Mechanisms detected | Mechanisms localized | False calls on non-disruptive variants")
for row in design_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("%from: docs/audit/COMAVI_v11_evidence_resolution.csv; reference_outputs/COMAVI_evidence_ledger.csv complex_amount_assay")
A("Ten binding expectations rest on assays that measure how much complex forms rather than binding itself: "
  "immunoprecipitation, pull-down, two-hybrid, TR-FRET and multiplexed interaction assays (S1 Table). A change in the amount of complex can come "
  "from a change in binding or in the stability of the assembly, so the six of these that expect a change (%s) are met by a "
  "call in the measured direction on either assembled axis (Methods). The four that expect no change (%s) are scored as "
  "binding expectations only, because a normal amount of complex cannot show that the assembled complex is stable. The rule "
  "is fixed by the assay recorded for each expectation and applies equally to COMAVI and to every structure-blind strategy. "
  "At 2.5 kcal/mol it changes three grades from inconsistent to consistent, for %s, raising the grade total from %s to %s "
  "(%s to %s). At 1.0, 1.5, 2.0 kcal/mol and the prediction-interval limit it adds %d, %d, %d and %d grade points. Without "
  "the rule, %d of the %d mechanisms would be localized rather than %d, and COMAVI's margin over calling no variant "
  "disrupted, a strategy the rule cannot affect, would be %s rather than %s. In direction agreement the six are scored as one joint row, "
  "and %d of them are met."
  % (names(scope_names), names(neutral_amount), names(sorted(vn(s, v) for s, v in zip(erp[erp.threshold == "t25"].system, erp[erp.threshold == "t25"].variant))),
     half(tot25 - 3), half(tot25), r(Fraction(int((tot25 - 3) * 2), 2 * N), 3), r(SCORE25, 3), gain_by_t["t10"], gain_by_t["t15"],
     gain_by_t["t20"], gain_by_t["tSAP"], att0["k"], att0["n"], att["k"], r(Fraction(int((tot25 - 3 - 27) * 2), 2 * N), 3),
     r(Fraction(int((tot25 - 27) * 2), 2 * N), 3), asm_either[0]))
A("")
A("## Section 4. Half credit for calls on axes the literature never tested")
A("")
A("%from: reference_outputs/COMAVI_rubric_sensitivity.json")
A("The primary grading treats a call that recovers no expected effect as inconsistent even when every change it asserts lies "
  "on an axis the literature never tested. As a sensitivity analysis only, such a call was credited 0.5 if it contradicted "
  "nothing measured. One grade changes at 2.5 kcal/mol, %s, and the score rises from %s to %s (95%% interval %s). The margin "
  "over the strongest structure-blind strategy, re-selected under the same rule, becomes %s (%s), and %s of resamples give a "
  "margin of zero or less. %s is credited at every operating point, and %s also at 1.0, 1.5 and 2.0 kcal/mol."
  % (names(chg_by["2.5"]), r(SCORE25, 3), r(Fraction(int((tot25 + 0.5) * 2), 2 * N), 3), ci(*rsr["score_ci95"]),
     r(Fraction(int((tot25 + 0.5 - 27) * 2), 2 * N), 3), ci(*rsr["margin_ci95"]),
     "%s%%" % r(exact_p(rsr["margin_p_at_or_below_zero"], rs["draws"]) * 100, 0),
     names(chg_by["2.5"]), names([x for x in chg_by["1.0"] if x not in chg_by["2.5"]])))
assert all(chg_by[k] == chg_by["2.5"] for k in ("2.5", "2.9/2.9/3.5")) and \
    all(set(chg_by[k]) == set(chg_by["2.5"]) | {"KRAS G12D"} for k in ("1.0", "1.5", "2.0"))
A("")
A("## Section 5. Operating-point diagnostics")
A("")
A("%from: reference_outputs/COMAVI_operating_point_diagnostics.json")
A("Table D gives, at each operating point, the mean grade on the mechanisms, the mean grade on the non-disruptive variants "
  "and the pooled score. Raising the threshold trades the first for the second. Over %s paired resamples of whole systems, "
  "the mechanism arm fell between 1.0 kcal/mol and the prediction-interval limit in %s, and the non-disruptive "
  "arm rose in %s. The pooled score weights the two arms by their share of the benchmark, %d mechanisms and %d "
  "non-disruptive variants (a structural fraction of %s). Reweighting the arms, 2.5 kcal/mol remains the pooled optimum for "
  "any weight on the mechanism arm between %s and %s, an interval that contains the benchmark's own fraction and the "
  "unweighted mean of the arms. Its place at the top is therefore not an artifact of the class balance, although the "
  "overlapping intervals of Fig 4A mean that the benchmark does not establish it as the best operating point."
  % (thousands(pb["n_draws"]), pct_of(pb["p_recovery_fell_t10_to_tSAP"]), pct_of(pb["p_rejection_rose_t10_to_tSAP"]),
     N_MECH, N_ND, r(Fraction(N_MECH, N), 3), r(cr["t25_optimal_for_structural_weight_above"], 3),
     r(cr["t25_optimal_for_structural_weight_below"], 3)))
A("")
A("@tablecap Table D. Mechanism arm, non-disruptive arm and pooled score at each operating point, with 95% intervals from resampling whole systems.")
A("@tabhead Operating point | Mechanism arm (29) | Non-disruptive arm (27) | Pooled score")
for row in op_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 6. Structural-context comparators and the interface bonus")
A("")
A("%from: reference_outputs/isds_v1/ISDS_v1_component_binary_baselines.csv; reference_outputs/COMAVI_tier_construction.json")
A("Interface occupancy enters the tier directly, through the interface bonus, so the tier was compared on the %d-variant "
  "prioritization population (%d mechanisms, %d non-disruptive variants in %d systems) with interface status alone and with "
  "the tier recalculated without its interface bonus (Table E; S3 Fig). The bonus reaches %d of the %d mechanisms, all six "
  "with mixed fold and interface expectations and %s of the interface destabilizers; the six it does not reach are %s. "
  "Removing it demotes four variants from Tier 2 to Tier 3, %s, of which only CALM1 D96V is non-disruptive. With the bonus, "
  "%d of the %d Tier 1 variants carry a structural mechanism, against %d of %d in Tier 2, %d of %d in Tier 3 and %d of %d "
  "in Tier 4. Within-system permutation p values are exact (S1 Text)."
  % (len(pvv), int(y46.sum()), int((~y46).sum()), pvv.system.nunique(), reach["structural_with_bonus"], reach["structural_total"],
     reach["by_class"]["ppi_destab_mechanism"].replace("/", " of "), names(no_bonus), names(dem), *tf["Tier 1"], *tf["Tier 2"],
     *tf["Tier 3"], *tf["Tier 4"]))
A("")
A("@tablecap Table E. Structural-context screens on the 46-variant prioritization population.")
A("@tabhead Screen | Mechanisms recovered (sensitivity) | Non-disruptive variants rejected (specificity) | Fisher p | Within-system permutation p")
for row in cmp_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 7. Evidence states")
A("")
A("%from: reference_outputs/isds_v1/ISDS_v1_evidence_state_tables.json; reference_outputs/COMAVI_tier_construction.json specificity_gain")
A("Crossing the energy calls with the tier (Tier 1\u20132 against Tier 3\u20134) gives four descriptive evidence states "
  "(Table F). At 2.5 kcal/mol, requiring the tier to agree with an energy call leaves sensitivity unchanged (%d of %d) and "
  "raises specificity from %s (%d of %d) to %s (%d of %d): the three non-disruptive variants called by energy alone, %s, "
  "fall in Tier 3 or 4. Resampling whole systems puts a 95%% interval of %s on that gain, and %d%% of resamples show no "
  "gain. The states are descriptive; none is used to grade a mechanism."
  % (c_["tp"], c_["tp"] + c_["fn"], r(Fraction(a_["tn"], a_["tn"] + a_["fp"]), 3), a_["tn"], a_["tn"] + a_["fp"],
     r(Fraction(c_["tn"], c_["tn"] + c_["fp"]), 3), c_["tn"], c_["tn"] + c_["fp"], names(eo),
     ci(*sg["reference_t25_cluster_bootstrap_ci"], nd=2), sg["fraction_resamples_no_gain_pct_reportable"]))
A("")
A("@tablecap Table F. Evidence states on the 46-variant prioritization population (structural mechanisms / non-disruptive variants).")
A("@tabhead Operating point | Convergent | Energy only | Context only | Neither")
for row in state_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 8. Alternative priority-score formulas and score bands")
A("")
A("%from: reference_outputs/COMAVI_priority_sensitivity_summary.json; reference_outputs/isds_v1")
A("The reported priority score reaches an AUC of %s (%s of the %d mechanism\u2013non-disruptive pairs ordered correctly; "
  "95%% interval %s from resampling whole systems, %s draws). Leaving one system out gives %s to %s, and the score exceeds "
  "both of its components in every fold. Table G compares nine alternative formulas, all computed on the same population. "
  "Five order slightly more pairs than the reported score and one orders the same number; the largest excess is %s pairs "
  "(%s, AUC %s), small against the width of the interval on the reported AUC. The alternatives that keep the graded context "
  "term rank the variants almost as the reported score does (Spearman \u03c1 %s to %s). One of them, the cohort-standardized combination, is not a per-variant score: "
  "it z-scores each term against the variants present, so a variant's value changes when others are added or removed, and "
  "it cannot be computed for a new variant without a reference cohort. It is retained only as a comparator. The reported "
  "score uses fixed constants and depends on no other variant, which is what makes it usable prospectively. Because the "
  "same benchmark informed the score's development, all of these are in-sample estimates."
  % (r(Fraction(int(rep_pairs * 2), 2 * npairs), 3), num(rep_pairs), npairs, ci(*isb["isds_roc_auc_ci95"], full=True),
     thousands(isb["draws_requested"]), r(min(loso_exact.values()), 3), r(max(loso_exact.values()), 3), num(top_alt[1] - rep_pairs),
     top_alt[0], r(Fraction(int(top_alt[1] * 2), 2 * npairs), 3), r(min(rho_ctx), 3), r(max(rho_ctx), 3)))
A("")
A("@tablecap Table G. Alternative priority-score formulas on the 46-variant prioritization population.")
A("@tabhead Formula | Graded context term | ROC AUC | Spearman \u03c1 with the reported score")
for row in alt_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("Partitioning the population into five equal score bands shows where the mechanisms sit (Table H): %d of the %d variants "
  "in the highest band carry a structural mechanism, against none of the %d in the lowest."
  % (int(bands.iloc[-1].n_structural_mechanism), int(bands.iloc[-1].n_variants), int(bands.iloc[0].n_variants)))
A("")
A("@tablecap Table H. Structural mechanisms by priority-score band.")
A("@tabhead Score band | Variants | Structural mechanisms | Fraction")
for row in band_rows:
    A("@tabrow " + " | ".join(row))
A("")
SOURCES["S2_Text"] = "\n".join(S2) + "\n"


# ================================================================ S3 Text
# S3 Text: pathogenicity comparison, robustness and independent validation.
# Executed after v11_si_text_edits.py, v11_si_s1.py and v11_si_s2.py in the same namespace.

# ---- Section 1: silent results by clinical label (benchmark population)
de = J("COMAVI_directional_enrichment.json")
t2 = de["table_2x2"]
assert de["cohort"]["n"] == N and de["labelled"]["n"] == 39 and de["threshold"] == 2.5
assert t2["firing_benign"] == 0 and sum(t2.values()) == de["labelled"]["n"]
_, p_de = fisher_exact([[t2["firing_pathogenic"], t2["firing_benign"]], [t2["silent_pathogenic"], t2["silent_benign"]]])
assert abs(p_de - de["fisher"]["p_value"]) < 5e-6
pbs = de["pathogenic_but_silent"]
assert pbs["n"] == t2["silent_pathogenic"] == len(pbs["variants"])
phen = {vn(s, v): str(p_) for s, v, p_ in zip(B.system, B.variant, B.phenotype)}
gof_silent = sorted(v for v in pbs["variants"] if phen[v] == "pathogenic_gof")
gof_all = sorted(v for v, p_ in phen.items() if p_ == "pathogenic_gof")
gof_firing = sorted(v for v in gof_all if v not in pbs["variants"])
assert set(gof_silent) | set(gof_firing) == set(gof_all)

# ---- Section 2: stress tests
sys_tot = B.assign(g=grades("t25")).groupby("system").g.agg(["sum", "count"])
loso_sc = {s: Fraction(int((tot25 - sys_tot.loc[s, "sum"]) * 2), 2 * int(N - sys_tot.loc[s, "count"])) for s in sys_tot.index}
lo_sys, hi_sys = min(loso_sc, key=loso_sc.get), max(loso_sc, key=loso_sc.get)
lcsv = pd.read_csv(RO / "stress_tests/comavi_leave_one_system_out.csv")
assert len(lcsv) == len(loso_sc) == 14
# n_dropped counts every scored row of the system, excluded candidates included;
# mc_n is the number of graded variants that remain, which is what the score is over.
for s_, n_, v_ in zip(lcsv.dropped, lcsv.mc_n, lcsv.mc):
    assert abs(float(loso_sc[s_]) - float(v_)) < 6e-5 and int(N - sys_tot.loc[s_, "count"]) == n_, s_
vb = J("stress_tests/COMAVI_stress_verified_statistics.json")
cl_mc, vb_mc = vb["system_cluster_ci95"]["mc"]["full_precision"], vb["variant_bootstrap_ci95"]["mc"]["full_precision"]
noise_mc, flips = z["noise_mc"], z["noise_flips"]
assert len(noise_mc) == 500
frag = pd.read_csv(RO / "stress_tests/comavi_noise_fragility.csv")
fragB = frag[[k in bkeys for k in zip(frag.system, frag.variant)]]
assert len(fragB) == N
unstable = fragB[fragB.flip_rate > 0].sort_values("flip_rate", ascending=False)
unstable_txt = ["%s (%s)" % (vn(s, v), r(Fraction(int(round(fr * 500)), 500), 3)) for s, v, fr in
                zip(unstable.system, unstable.variant, unstable.flip_rate)]
for fr in unstable.flip_rate:
    assert abs(fr * 500 - round(fr * 500)) < 1e-6
keep("noise", {"mean": float(noise_mc.mean()), "sd": float(noise_mc.std()), "flips_mean": float(flips.mean()),
               "flips_max": int(flips.max()), "unstable_calls": unstable_txt},
     "reference_outputs/stress_tests/comavi_stress_draws.npz; comavi_noise_fragility.csv (benchmark rows)")

# ---- Section 3: SKEMPI 2.0
skv = J("COMAVI_skempi_validation.json")
cx = skv["complexes"]
def pfmt(p_):
    """p to three decimals, or to two significant figures below 0.01."""
    if p_ >= 0.01:
        return r(p_, 3)
    nd = 1
    while p_ * 10 ** nd < 10:
        nd += 1
    return r(p_, nd)


sk_rows = []
for key, name in (("barnase_barstar_all", "Barnase\u2013barstar, all"),
                  ("barnase_barstar_excluding_glu73", "Barnase\u2013barstar, without Glu73"),
                  ("tem1_blip", "TEM1\u2013BLIP")):
    d = cx[key]
    sk_rows.append([name, "%d" % d["n"], "%s (%s)" % (r(d["pearson_r"], 2), pfmt(d["pearson_p"])),
                    "%s (%s)" % (r(d["spearman_rho"], 2), pfmt(d["spearman_p"])),
                    "%d of %d" % (d["sign_correct"], d["destabilizing_n"])])
for key, name in (("combined_including_glu73", "Both complexes"), ("combined_excluding_glu73", "Both complexes, without Glu73")):
    d = skv[key]
    sk_rows.append([name, "%d" % d["n"], "%s (%s)" % (r(d["pearson_r"], 2), pfmt(d["pearson_p"])),
                    "%s (%s)" % (r(d["spearman_rho"], 2), pfmt(d["spearman_p"])),
                    "%d of %d" % (d["sign_correct"], d["destabilizing_n"])])
assert cx["barnase_barstar_all"]["n"] + cx["tem1_blip"]["n"] == skv["combined_including_glu73"]["n"] == 48
g73 = skv["glu73_cluster"]
assert g73["all_measured_destabilizing"] and g73["all_predicted_stabilizing"] and g73["n"] == 7
bbv = pd.read_csv(RD / "supplement/skempi/bb_binding_validation.csv")
methods = bbv.method.value_counts().to_dict()

# ---- Section 4: rank x grade
rgo = J("isds_v1/ISDS_v1_rank_grade_overlap.json")
rec_ = rgo["recovery"]
assert (rgo["population"]["n_prioritized"], rgo["population"]["n_structural_truth"], rgo["population"]["top_k"]) == (46, 21, 20)
gin, gout = rgo["grade_breakdown_in_top_k"], rgo["grade_breakdown_outside_top_k"]
assert (gin["consistent"], gin["partial"], gin["inconsistent"]) == (8, 4, 5) and gout == {"consistent": 0, "partial": 0, "inconsistent": 4}
inc_in = {d["variant"]: d for d in rgo["inconsistent_inside_top_k"]}
assert set(inc_in) == {"R755W", "Q542L", "E545K", "N102T", "E542K"}
callB = {v: c for v, c in zip(B.variant, B.mech_t25)}
assert callB["R755W"] == callB["E545K"] == "PPI stabilization"
assert callB["E542K"] == "Multimer fold destabilization at interface"
assert callB["Q542L"] == callB["N102T"] == "Interface variant (DDG neutral)"
outs = sorted(rgo["structural_outside_top_k_detail"], key=lambda d: d["rank"])
assoc = rgo["association"]
r755 = B[B.variant == "R755W"].iloc[0]
cx_cols = [c for c in B.columns if c.startswith("ddg_fold_") and not c.endswith(("_sd", "_runs")) and "_ci95" not in c
           and "_distinguishable" not in c and "_indistinguishable" not in c and "_confident" not in c]
r755_cx = [float(r755[c]) for c in cx_cols if pd.notna(r755[c])]
assert len(r755_cx) == 1 and abs(r755_cx[0] - (-0.94)) < 0.005 and abs(r755.ddg_monomer - 0.97) < 0.005

# ---- Section 5: curation audit and correction history (wording checked against the records)
lb_gate = (RD / "scripts/verify_ledger_token_basis.py").read_text()
assert "def main" in lb_gate
wf = "\n".join(p_.read_text() for p_ in (RD / ".github/workflows").glob("*.yml"))
assert "verify_ledger_token_basis.py" in wf
v80 = (RD / "scripts/apply_ledger_corrections_v80.py").read_text()
for needle in ("WT 56 +/- 6, G12D 270 +/- 46, G12V 411 +/- 40", "REJECTED, evidence insufficient",
               "98 - 1 (G674R) + 2 (H718Y monomer, I500V binding) = 99", "Ollila 2006", "Hinrichsen 2013", "Lindsay 2025",
               "performs no co-immunoprecipitation"):
    assert needle in v80, needle
v79 = (RD / "scripts/apply_ledger_corrections_v79.py").read_text()
for needle in ("C697F", "G12D", "G12V", "13.6"):
    assert needle in v79, needle
assert round(270 / 56, 1) == 4.8 and round(411 / 56, 1) == 7.3
c697 = B[B.variant == "C697F"].iloc[0]
assert c697.expected_mech_class == SIL and c697.mech_consistency_t25 == "inconsistent" and c697.monomer_burial == "buried_core"
top_band = pvv[pvv.isds_v1 >= 0.8]
assert list(top_band[~top_band.structural_ground_truth.astype(bool)].variant) == ["C697F"]
assert pvv.sort_values("isds_v1", ascending=False).iloc[0].variant == "C697F"
kr = {v: B[B.variant == v].iloc[0] for v in ("G12D", "G12V")}
assert all(kr[v].mech_consistency_t25 == "inconsistent" for v in kr)
assert sorted(d["variant"] for d in outs) == ["G12D", "G12V", "H718Y", "I500V"]
assert all(B[B.variant == v].iloc[0].mech_t25 == "No structural effect detected" for v in ("G12D", "G12V", "H718Y", "I500V"))
mm = {d["variant"]: d for d in msd["variants"]}
meas_kras = {v: float(mm[v]["measured_kcal_mol"]) for v in ("G12D", "G12V")}
assert all(abs(mm[v]["predicted_binding"]) < 2.5 for v in ("G12D", "G12V"))
assert abs(meas_kras["G12D"] - 0.93) < 0.005 and abs(meas_kras["G12V"] - 1.18) < 0.005
g674 = can[can.variant == "G674R"].iloc[0]
# Before v8.0 (parent of commit c6f5ff6) G674R was curated structurally_silent and graded
# inconsistent on this same call; the call is unchanged, so the withdrawal removed a false call.
assert g674.mech_t25 == "Both fold destabilization" and g674.expected_mech_class == "structurally_uncommitted"
assert (g674.system, g674.variant) not in bkeys

# ---------------------------------------------------------------- source text
S3 = []
A = S3.append
A("@title S3 Text. Pathogenicity comparison, robustness and independent validation.")
A("")
A("@contents")
A("")
A("## Section 1. Silent results by clinical label")
A("")
A("%from: reference_outputs/COMAVI_directional_enrichment.json")
A("At the 2.5 kcal/mol reference, %d of the %d benchmark variants have at least one energy output that reaches the "
  "threshold and %d are silent; an output counts when a confidence-passing axis reaches the threshold, and binding values "
  "indistinguishable from zero are excluded. Of the %d variants with a pathogenic or benign classification (gain of function "
  "counted as pathogenic), all %d with an output at the threshold are pathogenic, whereas the %d silent variants include %d pathogenic and %d benign ones (Fisher "
  "exact p = %s). The %d pathogenic but silent variants span %d of the %d systems: %s. %s of them are curated as gain of "
  "function (%s), whose pathogenic effect is a change in activity rather than a loss of structure. Because the benchmark was selected for mechanism "
  "evidence, these fractions are not clinical predictive values; they establish only that the two kinds of result differ. "
  "No benign variant here had an output at the threshold, whereas silence occurred in both classes, so a silent result "
  "constrains nothing about pathogenicity."
  % (de["cohort"]["firing"], N, de["cohort"]["silent"], de["labelled"]["n"], t2["firing_pathogenic"],
     t2["silent_pathogenic"] + t2["silent_benign"], t2["silent_pathogenic"], t2["silent_benign"], r(p_de, 4),
     pbs["n"], pbs["systems_spanned"], pbs["systems_total"], names(sorted(pbs["variants"])), Cap(WORDS[len(gof_silent)]),
     names(gof_silent)))
A("")
A("## Section 2. Stress tests")
A("")
A("%from: verification/stress_tests.py; reference_outputs/stress_tests")
A("The main text reports the mechanism-pattern score with a 95%% interval from resampling whole systems (%s), because "
  "variants within a system share a structure and a curation pass. Resampling variants instead, which ignores that "
  "clustering, gives %s. Leaving one system out at a time gives values from %s to %s (Table A). Adding Gaussian noise at the standard error of each five-replicate FoldX mean, %s times, gave a mean score of %s "
  "(standard deviation %s); a draw changed %s calls on average and at most %d. Of the %d benchmark variants, %d kept the same "
  "call in every draw; the others, with the fraction of draws in which their call changed, are %s."
  % (ci(*cl_mc, full=True), ci(*vb_mc, full=True), r(loso_sc[lo_sys], 3), r(loso_sc[hi_sys], 3), thousands(len(noise_mc)), r(float(noise_mc.mean()), 3, full=True), r(float(noise_mc.std()), 3, full=True),
     r(float(flips.mean()), 2, full=True), int(flips.max()), N, N - len(unstable), names(unstable_txt)))
A("")
A("@tablecap Table A. Mechanism-pattern score with each system left out, at 2.5 kcal/mol.")
A("@tabhead System left out | Variants removed | Score on the remaining variants")
for s_ in sorted(loso_sc, key=lambda s_: (float(loso_sc[s_]), SYS[s_])):
    A("@tabrow %s | %d | %s (%s of %d)" % (SYS[s_], int(sys_tot.loc[s_, "count"]), r(loso_sc[s_], 3),
                                          half(tot25 - sys_tot.loc[s_, "sum"]), int(N - sys_tot.loc[s_, "count"])))
A("")
A("## Section 3. Independent binding-axis validation on SKEMPI 2.0")
A("")
A("%from: reference_outputs/COMAVI_skempi_validation.json; supplement/skempi")
A("The hemoglobin series is a high-quality within-system anchor but does not test unrelated interfaces. Two complexes from "
  "SKEMPI 2.0 [@jankauskaite2019] were therefore scored with the frozen workflow: barnase\u2013barstar (%d substitutions) "
  "and TEM1\u2013BLIP (%d substitutions on the BLIP side, where residue numbering maps unambiguously). "
  "Measured binding \u0394\u0394G is RT ln(K_D,mut/K_D,wt). Table B gives the agreement for each complex. Seven barnase "
  "substitutions at Glu73 (measured %s to %s kcal/mol, all destabilizing) are all predicted to stabilize binding (%s to %s "
  "kcal/mol), a single-site sign inversion. They are kept in the primary statistics, and the fit without them, which was "
  "introduced with the analysis rather than prespecified, is reported as a sensitivity view."
  % (cx["barnase_barstar_all"]["n"], cx["tem1_blip"]["n"], r(g73["measured_min"], 2), r(g73["measured_max"], 2), r(g73["predicted_min"], 2), r(g73["predicted_max"], 2)))
A("")
A("@tablecap Table B. FoldX binding \u0394\u0394G against SKEMPI 2.0 measurements.")
A("@tabhead Complex | Substitutions | Pearson r (p) | Spearman \u03c1 (p) | Sign correct on measured destabilizations")
for row in sk_rows:
    A("@tabrow " + " | ".join(row))
A("")
A("## Section 4. Do high-priority variants also have their mechanism graded correctly?")
A("")
A("%from: reference_outputs/isds_v1/ISDS_v1_rank_grade_overlap.json")
A("The priority score and the mechanism grade answer different questions, and a variant ranked highly but graded "
  "inconsistent would point a user at the right variant with the wrong mechanistic hypothesis. On the %d-variant "
  "prioritization population, %d of the %d structural mechanisms are ranked in the top %d: %d are graded consistent, %d "
  "partial and %d inconsistent. Of the five inconsistent, two are placed on the wrong axis, MLH1 R755W (rank %d), a fold "
  "mechanism called as strengthened binding, and PIK3CA E542K (rank %d), an interface destabilizer called as instability "
  "of the assembly; one, PIK3CA E545K (rank %d), is called on the expected binding axis but in the opposite direction; and "
  "two, MLH1 Q542L (rank %d) and HBB N102T (rank %d), receive no energy call and are ranked by their structural context. "
  "For R755W the literature reports reduced abundance with binding preserved, while neither stability prediction reaches "
  "the threshold (%s kcal/mol for the isolated subunit and %s in the complex). R755W and E542K would therefore send a user to "
  "the right protein with the wrong experiment. The four mechanisms outside the top %d, %s, are all graded inconsistent, so "
  "no fully resolved mechanism is ranked below the cut. Across the %d mechanisms, the priority score and the grade are "
  "correlated (Spearman \u03c1 = %s, p = %s), and the mean grade is %s inside the top %d against %s outside it (one-sided "
  "Mann\u2013Whitney p = %s). With only four mechanisms outside the cut, the contrast is descriptive."
  % (rgo["population"]["n_prioritized"], rec_["structural_in_top_k"], rgo["population"]["n_structural_truth"], rgo["population"]["top_k"],
     gin["consistent"], gin["partial"], gin["inconsistent"], inc_in["R755W"]["rank"], inc_in["E542K"]["rank"], inc_in["E545K"]["rank"],
     inc_in["Q542L"]["rank"], inc_in["N102T"]["rank"], r(r755.ddg_monomer, 2), r(r755_cx[0], 2), rgo["population"]["top_k"],
     names(["%s %s (rank %d)" % (d["gene"].upper(), d["variant"], d["rank"]) for d in outs]), rgo["population"]["n_structural_truth"],
     r(assoc["spearman_rho_score_vs_grade"], 3), r(assoc["spearman_p"], 4), r(rgo["mean_grade"]["in_top_k"], 3),
     rgo["population"]["top_k"], r(rgo["mean_grade"]["outside_top_k"], 3), r(assoc["mannwhitney_p_one_sided"], 3)))
A("")
A("## Section 5. Curation audit and correction history")
A("")
A("%from: scripts/verify_ledger_token_basis.py; scripts/apply_ledger_corrections_v79.py; scripts/apply_ledger_corrections_v80.py")
A("Every committed expectation carries a one-line evidence basis and a primary citation beside its expected direction (S1 "
  "Table). A gate run in continuous integration, scripts/verify_ledger_token_basis.py, reads each basis, keeps the clauses "
  "that use the vocabulary of the axis being committed, and fails if the direction asserts disruption where those clauses "
  "assert that the axis is intact, or the reverse.")
A("")
A("When the gate was introduced, three of the then 98 committed expectations failed it. All three were binding expectations "
  "that one earlier curation review had changed from neutral to destabilizing without changing the bases beside them: MSH2 "
  "C697F, whose basis read \u201cATPase catalytic defect, protein stable, interaction intact\u201d, and KRAS G12D and G12V, "
  "whose bases recorded wild-type-like RAF1 binding with impaired GTP hydrolysis and noted that residue 12 lies 13.6 \u00c5 "
  "from the RAF1 interface. The three were returned to neutral, the direction their bases stated.")
A("")
A("The next revision read the cited sources directly and found that the KRAS bases themselves misread their source: Hunter "
  "et al. [@hunter2015] report RAF affinities of 270 nM for G12D and 411 nM for G12V against 56 nM for the wild type, a "
  "4.8-fold and a 7.3-fold loss. Basis and direction were corrected together, and both are again destabilizing binding "
  "expectations. A proposal to make MSH2 C697F destabilizing was rejected: its source could not be read, and the neutral "
  "commitment is supported by a source that could, a co-immunoprecipitation study that found the interaction intact "
  "[@ollila2006]. The same revision withdrew the binding expectation of MSH2 G674R, whose basis cited a study that performs "
  "no co-immunoprecipitation (the variant now has no committed expectation and is not part of the benchmark), and committed "
  "two expectations from the sources' own measurements: isolated-subunit destabilization for MLH1 H718Y [@hinrichsen2013] "
  "and strengthened SMAD3 binding for SMAD4 I500V [@lindsay2025]. The ledger then held 99 expectations, of which the 97 "
  "belonging to benchmark variants are in S1 Table. Each change, with its reason, is recorded in the correction scripts of "
  "the repository.")
A("")
A("These corrections are reported rather than absorbed because their effects run in both directions. Withdrawing the "
  "G674R expectation removed from the benchmark a variant that had been curated as non-disruptive and on which COMAVI "
  "predicts fold destabilization (%s kcal/mol for the isolated subunit), so it removed a false call. The others lower the "
  "score. KRAS G12D and G12V are "
  "graded inconsistent, because COMAVI predicts no binding change at the reference and their measured losses (%s and %s "
  "kcal/mol) are themselves below it. MLH1 H718Y and SMAD4 I500V, committed from their sources' measurements, are mechanisms "
  "that no axis detects. These four are the mechanisms ranked outside the top 20 (Section 4). MSH2 C697F remains a curated "
  "negative on which COMAVI predicts %s kcal/mol of isolated-subunit destabilization at a buried-core residue, so it is a "
  "false call and the highest-ranked variant in the prioritization population; the rejected change would have raised the "
  "score. A gate of "
  "this kind cannot detect every curation error. It catches contradictions that are written down, but not a basis that "
  "misreads the source it names, is silent on the axis, cites the wrong article, or records a mistaken reading faithfully. "
  "A passing gate means that no committed direction contradicts its own recorded evidence, which is weaker than every "
  "direction being correct." % (r(g674.ddg_monomer, 2), r(meas_kras["G12D"], 2), r(meas_kras["G12V"], 2), r(c697.ddg_monomer, 1)))
A("")
A("@references")
SOURCES["S3_Text"] = "\n".join(S3) + "\n"
for k_, v_ in SOURCES.items():
    check_keys(v_)


# ---------------------------------------------------------------- write / check
AUDIT = RD / "docs" / "audit" / "COMAVI_v11_si_text_values.json"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="regenerate in memory and compare with the committed sources; write nothing")
    args = ap.parse_args()
    audit = json.dumps(VALUES, ensure_ascii=False, indent=1, sort_keys=True, default=str) + "\n"
    targets = {SI_SRC / ("%s_source.txt" % k): v for k, v in SOURCES.items()}
    targets[AUDIT] = audit
    if args.check:
        bad = [p.relative_to(RD) for p, v in targets.items() if not p.exists() or p.read_text() != v]
        for b in bad:
            print("FAIL: %s differs from a fresh generation" % b)
        if bad:
            return 1
        print("PASS: S1-S3 Text sources and the values audit regenerate from the records (%d values)" % len(VALUES))
        return 0
    SI_SRC.mkdir(parents=True, exist_ok=True)
    for p, v in targets.items():
        p.write_text(v)
        print("wrote %s (%d characters)" % (p.relative_to(RD), len(v)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
