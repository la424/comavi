# PROVENANCE RECORD, NOT A RUNNABLE STEP OF THE PIPELINE.
# This script was applied once, to the v10.1 manuscript source (commit e1dcbb3),
# to write the v11 rewrite; the values it read are in
# docs/audit/COMAVI_v11_manuscript_values.json. Every later change to the text
# is listed, with its reason, in docs/audit/COMAVI_v11_post_draft_text_edits.csv.
# Each edit locates its target by a unique prefix of a v10.1 line and asserts it,
# so re-running it against a later source stops without writing anything.
#!/usr/bin/env python3
"""v11 manuscript edits, with every number computed from the committed records.

    python v11_manuscript_edits.py --repo /path/to/comavi [--dry-run]

Each edit replaces one whole source line (found by a unique prefix), deletes
one paragraph, inserts paragraphs after an anchor, or replaces one exact
substring that must occur exactly once. Every value is read from a committed
record at run time, and every premise a sentence states (an ordering, a
"half", "at every threshold") is asserted first, so the script refuses to write
a sentence the records do not support. The values used are written to
docs/audit/COMAVI_v11_manuscript_values.json.
"""
import argparse
import importlib.util
import json
import pathlib
import re
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--repo", required=True)
ap.add_argument("--dry-run", action="store_true")
ARGS = ap.parse_args()
REPO = pathlib.Path(ARGS.repo)
RO = REPO / "reference_outputs"
SRC = REPO / "docs" / "manuscript" / "COMAVI_manuscript_source.txt"
sys.path[:0] = [str(REPO / "scripts"), str(REPO)]
import benchmark_population as bp  # noqa: E402
import apply_concordance_v5 as ac  # noqa: E402

_spec = importlib.util.spec_from_file_location("bda", REPO / "scripts" / "build_detection_vs_attribution.py")
bda = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bda)

J = lambda p: json.loads((RO / p).read_text())
TAGS = ["t10", "t15", "t20", "t25", "tSAP"]
TLAB = {"t10": "1.0", "t15": "1.5", "t20": "2.0", "t25": "2.5", "tSAP": "the prediction-interval limit"}
GRADE = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
SIL = "structurally_silent"
MINUS = "\u2212"
# Rounding rule for every reported value: strip float noise (10 decimals), then
# round half to even (ISO 80000-1). Formatting a float directly with "%.3f"
# resolves exact ties by binary representation instead of by a rule (v11: the
# weak-evidence margin bound is exactly -7/80, carried as -0.08749999999999997).
import decimal as _dec


def _fmt(x, places):
    q = _dec.Decimal(repr(round(float(x), 10))).quantize(_dec.Decimal(1).scaleb(-places), rounding=_dec.ROUND_HALF_EVEN)
    return ("%s" % q).replace("-", MINUS)


f3 = lambda x: _fmt(x, 3)
f2 = lambda x: _fmt(x, 2)
assert (f3(0.4375), f3(0.5125), f3(-0.08749999999999997), f3(34 / 53), f2(0.1875)) == \
    ("0.438", "0.512", MINUS + "0.088", "0.642", "0.19")
W = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
     11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen", 17: "seventeen",
     18: "eighteen", 19: "nineteen", 20: "twenty"}
Wc = {k: v.capitalize() for k, v in W.items()}
V = {}

# ------------------------------------------------------------------ records
can = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
B = bp.benchmark(can).reset_index(drop=True)
sw = {r["threshold_tag"]: r for r in J("COMAVI_threshold_sweep.json")["rows"]}
dva = J("COMAVI_detection_vs_attribution.json")
rb = J("COMAVI_rubric_baselines.json")
# Evidence-sensitivity values at full precision from the generator itself; the
# committed record holds them at 4 decimals, and re-rounding those can move the
# third decimal. Every full-precision value must round to the record.
_es_spec = importlib.util.spec_from_file_location("build_evidence_sensitivity",
                                                  REPO / "scripts" / "build_evidence_sensitivity.py")
_es_mod = importlib.util.module_from_spec(_es_spec)
_es_spec.loader.exec_module(_es_mod)
es = _es_mod.build(digits=None)


def _same4(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_same4(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(_same4(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return round(a, 4) == b
    return a == b


assert _same4(es, J("COMAVI_evidence_sensitivity.json")), "full-precision build disagrees with the committed record"
stress = J("stress_tests/COMAVI_stress_verified_statistics.json")
st = pd.read_csv(RO / "stress_tests" / "comavi_stress_tests.csv")
ccc = J("COMAVI_cohort_composition_checks.json")
rgo = J("isds_v1/ISDS_v1_rank_grade_overlap.json")
pc = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_pathogenicity_comparison.csv").set_index("descriptor")
cs = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_common_set.csv")
pv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
led = pd.read_csv(RO / "COMAVI_evidence_ledger.csv")
_r, FR = bda.build(can, return_frame=True)


def one(df, **kw):
    m = pd.Series(True, index=df.index)
    for k, v in kw.items():
        m &= df[k] == v
    assert int(m.sum()) == 1, (kw, int(m.sum()))
    return df[m].iloc[0]


# ------------------------------------------------------------------ populations
sil = (B.expected_mech_class == SIL).to_numpy()
br = (B.system == "brca1_brct").to_numpy()
V.update(n=len(B), systems=B.system.nunique(), mech=int((~sil).sum()), nondis=int(sil.sum()),
         n_cx=int((~br).sum()), cx_systems=int(B[~br].system.nunique()),
         cx_mech=int((~sil & ~br).sum()), cx_nondis=int((sil & ~br).sum()),
         n_br=int(br.sum()), br_mech=int((~sil & br).sum()), br_nondis=int((sil & br).sum()))
assert (V["n"], V["systems"], V["n_cx"], V["cx_systems"], V["n_br"]) == (bp.N_BENCHMARK, 14, len(pv), 13, 10), V
assert set(B[~sil & br].expected_mech_class) == {"fold_mechanism"}
excluded = sorted(set(zip(can.system, can.variant)) - set(zip(B.system, B.variant)))
assert sorted(v for _, v in excluded) == sorted(["R162W", "W117R", "G674R", "R1699L", "R1699Q"]), excluded
nh = int((B.gene.astype(str).str.upper() == "HBB").sum())
assert nh == 6 and {"I62V", "A1381T"} <= set(B.variant)
V["uniprot_named"] = V["n"] - nh - 2

# ------------------------------------------------------------------ pathogenicity comparison
lab = B.set_index(["system", "variant"])
cs["struct"] = [lab.loc[(s, v), "expected_mech_class"] != SIL for s, v in zip(cs.system, cs.variant)]
P = cs.clinical_y == 1
V.update(n_am=len(cs), P=int(P.sum()), Bn=int((~P).sum()), P_nondis=int((P & ~cs.struct).sum()),
         P_struct=int((P & cs.struct).sum()), B_nondis=int((~P & ~cs.struct).sum()))
assert sorted(cs[~P & cs.struct].variant) == ["H718Y"]
assert 2 * V["P_nondis"] == V["P"], "the text says half of the pathogenic class"
am_keys = set(zip(cs.system, cs.variant))
assert am_keys <= set(zip(B[~br].system, B[~br].variant)), "must nest in the prioritization population"
nc = B[(~br) & pd.Series([(s, v) not in am_keys for s, v in zip(B.system, B.variant)]).to_numpy()]
V["not_compared"] = len(nc)
assert V["n_am"] + V["not_compared"] == V["n_cx"] and set(nc.phenotype) <= {"vus", "unclassified"}
AUC = {k: float(pc.loc[d, "pathogenicity_auc"]) for k, d in
       [("am", "AM pathogenicity"), ("pr", "isds_v1"), ("en", "isds_energy_component"), ("cx", "isds_context_component")]}
assert int(pc.loc["AM pathogenicity", "n"]) == V["n_am"]
assert AUC["am"] > AUC["en"] > AUC["pr"] > AUC["cx"], AUC
pm_isds = pd.read_csv(RO / "isds_v1" / "ISDS_v1_primary_metrics.csv").set_index("score").loc["isds_v1"]
assert f3(pm_isds["roc_auc"]) == "0.834" and int(pm_isds["n"]) == 46, "structural AUC quoted in Results and Discussion"
V["auc"] = AUC
AMB = (0.34, 0.564)
amr = {v: one(cs, variant=v) for v in ("E6V", "R78G", "H1047R", "N127S")}
assert [int(amr[v].clinical_y) for v in ("E6V", "R78G", "H1047R", "N127S")] == [1, 1, 1, 0]
assert all(AMB[0] <= amr[v]["AM pathogenicity"] < AMB[1] for v in ("R78G", "H1047R", "N127S"))
assert sorted(cs[(cs.clinical_y == 0) & (cs["AM pathogenicity"] >= AMB[0]) & (cs["AM pathogenicity"] < AMB[1])].variant) == ["N127S"]
h = one(B, variant="H1047R")
assert "2" in str(h.comavi_tier) and not bool(FR.set_index(["system", "variant"]).loc[(h.system, "H1047R"), "any_fires"])
V["am_vals"] = {v: (float(amr[v]["AM pathogenicity"]), float(amr[v]["isds_v1"])) for v in amr}

# ------------------------------------------------------------------ evidence ledger
keys = set(zip(B.system, B.variant))
L = led[[(s, v) in keys for s, v in zip(led.system, led.variant)]].copy()
L["cls"] = L.evidence_type.str[:2]
tok = L.expected_token.astype(str).str.lower()
silk = set(zip(B[sil].system, B[sil].variant))
on_nd = pd.Series([(s, v) in silk for s, v in zip(L.system, L.variant)], index=L.index)
cc = L.cls.value_counts().to_dict()
weak = L.cls.isin(["E4", "E5"])
V.update(ledger_n=len(L), e1=cc["E1"], e2=cc["E2"], e3=cc["E3"], e4=cc["E4"], e5=cc["E5"],
         destab=int((tok == "destab").sum()), destab_quant=int(((tok == "destab") & L.cls.isin(["E1", "E2"])).sum()),
         e45=int(weak.sum()), e45_nd=int((weak & on_nd).sum()))
assert V["e1"] + V["e2"] == 40, "the evidence paragraph opens with 'Forty'"
assert int((weak & (tok == "neutral")).sum()) == V["e45"], "all population and inferred expectations neutral"
nolead = sorted(v for s, v in keys if (s, v) not in set(zip(L.system, L.variant)))
assert nolead == sorted(["D96V", "N98S", "F142L", "I500T"]), nolead

# ------------------------------------------------------------------ detection
t25, t10, tS = sw["t25"], sw["t10"], sw["tSAP"]
fc = lambda r: r["correct_rejection_n"] - r["correct_rejection_k"]
fci = lambda r: r["correct_rejection_n"] - r["correct_rejection_isolated_k"]
V.update(det_iso=t25["detection_isolated_k"], det=t25["detection_k"], att=t25["attribution_k"],
         rej=t25["correct_rejection_k"], added_fc=fc(t25) - fci(t25),
         cx_iso=t25["detection_complex_isolated_k"], cx_all=t25["detection_complex_k"],
         fc10=fc(t10) - fci(t10), fcS=fc(tS) - fci(tS))
# Sentences of the abstract and author summary that are kept as written.
# v11 final (evidence resolution): attribution 15 -> 18.
assert (V["det_iso"], V["det"], V["att"], V["added_fc"], V["mech"], V["nondis"]) == (7, 20, 18, 3, 29, 27)
assert rgo["recovery"]["structural_in_top_k"] == 17 and rgo["population"]["n_structural_truth"] == 21
assert all(sw[t]["detection_k"] > sw[t]["attribution_k"] for t in TAGS), "detection exceeds localization at every threshold"
# "Localizes most of what it detects": the wrong-axis cases, the mechanisms the
# either-axis rule credits, and the interval on the wrong-axis fraction.
_wrong = FR[FR.expected_mech_class.isin(bda.STRUCTURAL) & FR.any_fires.astype(bool)
            & ~FR.right_axis.map(lambda v: v is True)]
V["wrong_axis"] = sorted(zip(_wrong.gene.astype(str).str.upper(), _wrong.variant))
assert V["wrong_axis"] == [("MLH1", "R755W"), ("PIK3CA", "E542K")], V["wrong_axis"]
_c0 = can.copy()
_c0["assembled_expectation_scope"] = ""
_r0, _F0 = bda.build(_c0, return_frame=True)
_ra = lambda f: f.set_index(["system", "variant"]).right_axis.map(lambda v: v is True)
_cred = _ra(FR) & ~_ra(_F0)
V["credited_by_rule"] = sorted(v for (s, v), x in _cred.items() if x)
assert V["credited_by_rule"] == ["D351H", "L749P", "V11G"], V["credited_by_rule"]
assert _r0["attribution"]["k"] + len(V["credited_by_rule"]) == V["att"]
_scope = can[can.assembled_expectation_scope.fillna("") == ac.EITHER_ASSEMBLED_AXIS]
V["n_scope"] = len(_scope)
assert V["n_scope"] == 6 and set(_scope.expected_mech_class) <= {"ppi_destab_mechanism", "ppi_stab_mechanism"}
V["loc_gap_ci"] = dva["system_cluster"]["localization_gap"]["ci95"]
# The missed mechanisms (no axis reaches 2.5 kcal/mol), with measured magnitudes
# from committed records only (COMAVI_detection_vs_attribution.json).
MS = dva["missed_at_reference"]
assert MS["n_missed"] == V["mech"] - V["det"] == 9
_mv = {x["variant"]: x for x in MS["variants"]}
_below = sorted(v for v, x in _mv.items() if x["measured_kcal_mol"] is not None and x["measured_kcal_mol"] < 2.5)
_large = sorted(v for v, x in _mv.items() if x["measured_kcal_mol"] is not None and x["measured_kcal_mol"] >= 2.5)
_none = sorted(v for v, x in _mv.items() if x["measured_kcal_mol"] is None)
assert _below == ["G12D", "G12V", "L1664P", "R1751Q"] and _large == ["N102T", "V1736A"] \
    and _none == ["H718Y", "I500V", "Q542L"], (_below, _large, _none)
for v in _large:   # large measured effect, every predicted axis below the reference
    assert all(p is None or abs(p) < 2.5 for p in (_mv[v]["predicted_monomer"], _mv[v]["predicted_complex_fold"],
                                                   _mv[v]["predicted_binding"])), v
V["below_rng"] = (min(_mv[v]["measured_kcal_mol"] for v in _below), max(_mv[v]["measured_kcal_mol"] for v in _below))
V["mv"] = _mv
# Stringency: detection falls, the wrong-axis count shrinks, rejection rises.
_det = [sw[t]["detection_k"] for t in TAGS]
_att = [sw[t]["attribution_k"] for t in TAGS]
_rej = [sw[t]["correct_rejection_k"] for t in TAGS]
_wa = [d - a for d, a in zip(_det, _att)]
assert all(x >= y for x, y in zip(_det, _det[1:])) and _det[0] > _det[-1], _det
assert all(x >= y for x, y in zip(_wa, _wa[1:])) and _wa[0] > _wa[-1], _wa
assert all(x < y for x, y in zip(_rej, _rej[1:])), _rej
assert (_det[0] - _det[-1]) > (_att[0] - _att[-1]), "detection falls more than attribution"
V.update(det_t=_det, att_t=_att, rej_t=_rej, wa_t=_wa)
# Per-axis thresholds (COMAVI_per_axis_threshold_test.json, which calls the same build()).
PA = J("COMAVI_per_axis_threshold_test.json")
assert PA["reference_uniform"]["attribution_k"] == V["att"]
assert PA["gain_over_reference"]["attribution"] == 0 and PA["gain_over_reference"]["interface_attribution"] == 0
_ni = PA["interface_not_attributed_at_reference"]
assert PA["unreachable_by_any_threshold"]["k"] == _ni["k"] == 5, "every unattributed interface mechanism is unreachable"
_lost = max(PA["reference_uniform"]["correct_rejection_k"] - r["correct_rejection_k"] for r in PA["binding_threshold_sweep"])
V.update(pa_if=(PA["reference_uniform"]["interface_attribution_k"], PA["reference_uniform"]["interface_attribution_n"]),
         pa_ni=_ni["k"], pa_med=PA["unreachable_by_any_threshold"]["median_predicted_binding"], pa_lost=_lost)
assert _lost == 2 and V["pa_if"] == (7, 12)
asm = {t: sw[t]["detection_k"] - sw[t]["detection_isolated_k"] for t in TAGS}
assert all(v > 0 for v in asm.values())
assert all(fc(sw[t]) >= fci(sw[t]) for t in TAGS)
g = FR[FR.expected_mech_class.isin(bda.STRUCTURAL)]
only = g[g.any_fires.astype(bool) & ~g.mono_fires.astype(bool)]
assert len(only) == V["det"] - V["det_iso"]
V["n_sys_added"] = int(only.system.nunique())
pcls = dva["per_class"]
assert pcls["ppi_destab_mechanism"]["fires_only_assembled_axes"] + pcls["mixed_structural"]["fires_only_assembled_axes"] \
    > pcls["fold_mechanism"]["fires_only_assembled_axes"], "added detections mainly interface and mixed"
ifd = FR[FR.expected_mech_class == "ppi_destab_mechanism"]
V["n_if_destab"] = len(ifd)
assert (len(ifd), int(ifd.bind_fires.sum()), int(ifd.fold_fires.sum())) == (11, 4, 7)
SC = dva["system_cluster"]
D = dva["axis_designs"]
ib, al = D["isolated_subunit_and_binding"], D["all_three"]
olist = D["detected_only_through_assembled_stability"]
n_only, n_only_if = len(olist), sum(d["expected_class"] == "ppi_destab_mechanism" for d in olist)
assert (n_only, n_only_if) == (5, 3), (n_only, n_only_if)
V.update(d_all=al["detection"], d_ib=ib["detection"], a_all=al["attribution"], a_ib=ib["attribution"],
         fc_saved=al["false_calls"] - ib["false_calls"])
assert V["fc_saved"] == 2
for v in ("R78G", "R1334Q"):
    assert one(B, variant=v).mech_consistency_t25 == "partial", v
# Of the mechanisms detected only through assembled-complex stability, which are
# placed on the expected axis (v11 final: all but PIK3CA E542K).
_ok = _ra(FR)
_only_att = [d["variant"] for d in olist if _ok.loc[(d["system"], d["variant"])]]
_only_wrong = sorted("%s %s" % (d["gene"], d["variant"]) for d in olist if not _ok.loc[(d["system"], d["variant"])])
assert len(_only_att) == 4 and _only_wrong == ["PIK3CA E542K"], (_only_att, _only_wrong)
V.update(only_att=len(_only_att))
# Interface destabilizers placed on the expected axis.
V["if_att"] = int(ifd.right_axis.map(lambda v: v is True).sum())
assert V["if_att"] == 7, V["if_att"]

# ------------------------------------------------------------------ mechanism-pattern score
ev = rb["every_constant_label"]
best = ev["best_score"]
V.update(k=t25["whole_variant_k"], score=t25["whole_variant"], n_labels=ev["n_labels_graded"], best=best,
         all3=rb["baselines"]["always_fold_and_ppi"]["score"], margin=ev["margin_over_best_constant"])
assert set(ev["best_labels"]) == {"No structural effect detected", "Interface variant (DDG neutral)",
                                  "Structural variant \u2014 burial-driven (DDG neutral)",
                                  "Structural variant \u2014 contact-driven (DDG neutral)"}, ev["best_labels"]
assert abs(best - V["nondis"] / V["n"]) < 1e-3, "the strongest strategy scores the non-disruptive fraction"
assert tuple(ac.MECH_LABEL_ASSERTIONS["Both fold + PPI destabilization"][:3]) == (True, True, True)
assert abs(ev["labels"]["Both fold + PPI destabilization"]["score"] - V["all3"]) < 1e-9
prow = st[st.test == "permutation null (MC)"].iloc[0]
V.update(perm_mean=float(prow.null_mean), perm_sd=float(prow.null_sd), perm_p=float(prow.p_value))
mt = {t: sw[t]["whole_variant"] - best for t in TAGS}
assert min(mt.values()) > 0, "margin positive at every threshold"
bysys = es["by_system"]
lm = [v["margin_without_system"] for v in bysys.values()]
assert min(lm) > 0, "margin positive whichever system is left out"
contrib = sorted(((v["contribution_vs_calling_nothing"], s) for s, v in bysys.items()), reverse=True)
total = sum(c for c, _ in contrib)
assert [s for _, s in contrib[:2]] == ["brca1_brct", "hemoglobin_tetramer"], contrib[:3]
# v11 final: the two largest contributions are most of the margin, not more than
# all of it (contributions are in grade points; their sum is the margin x n).
assert contrib[0][0] + contrib[1][0] < total and 2 * (contrib[0][0] + contrib[1][0]) > total, (contrib[:2], total)
assert abs(total - (V["k"] - V["nondis"])) < 1e-9, (total, V["k"], V["nondis"])   # margin in grade points
V.update(c_brct=contrib[0][0], c_hb=contrib[1][0], c_total=total)
n_neg = sum(c < 0 for c, _ in contrib)
FB = es["views"]["full_benchmark"]
assert abs(FB["margin"] - V["margin"]) < 1e-3
V.update(lm=(min(lm), max(lm)), m_ci=FB["margin_ci95"], m_p=FB["margin_p_at_or_below_zero"], n_neg=n_neg)
# Leave-one-system-out range, exact from the graded benchmark (the stress CSV
# stores 4 decimals: 34/53 = 0.641509 is stored as 0.6415). Cross-checked
# against that CSV.
_g25 = B.mech_consistency_t25.map(GRADE)
_loso = [float(_g25[B.system != s_].mean()) for s_ in B.system.unique()]
lo_rng = [min(_loso), max(_loso)]
_lo_csv = [float(x) for x in st[st.test == "leave-one-system-out (MC range)"].iloc[0].null_mean.split("-")]
assert all(abs(a - b) < 5e-5 for a, b in zip(lo_rng, _lo_csv)), (lo_rng, _lo_csv)
V.update(ci=stress["system_cluster_ci95"]["mc"]["full_precision"], loso=lo_rng)

arms = {}
for t in TAGS:
    gg = B["mech_consistency_%s" % t].map(GRADE)
    assert gg.notna().all(), t
    arms[t] = (float(gg[~sil].mean()), float(gg[sil].mean()))
assert abs(arms["t25"][0] - rb["observed"]["arm_committed_structural"]) < 1e-3
assert abs(arms["t25"][1] - rb["observed"]["arm_curated_silent"]) < 1e-3
attr = {t: sw[t]["attribution_k"] / sw[t]["attribution_n"] for t in TAGS}
gaps = {t: attr[t] - arms[t][0] for t in ("t10", "t15", "t20", "t25")}
assert all(x > 0 for x in gaps.values()), ("structural arm below attribution at every uniform threshold", gaps)
unw = {t: (arms[t][0] + arms[t][1]) / 2 for t in TAGS}
assert max(TAGS, key=lambda t: sw[t]["whole_variant"]) == "t25" and max(TAGS, key=unw.get) == "t25", unw
d_arm = [arms[t][0] - arms[t][1] for t in TAGS]
cross = [i for i in range(4) if d_arm[i] > 0 >= d_arm[i + 1]]
assert len(cross) == 1 and all(x > 0 for x in d_arm[:cross[0] + 1]) and all(x <= 0 for x in d_arm[cross[0] + 1:]), d_arm
a_, b_ = TAGS[cross[0]], TAGS[cross[0] + 1]
V["cross"] = ("%s and %s kcal/mol" % (TLAB[a_], TLAB[b_]) if b_ != "tSAP"
              else "2.5 kcal/mol and the prediction-interval limit")
V.update(arms=arms, gaps=gaps, unweighted=unw)
g25 = B.mech_consistency_t25.map(GRADE)
V.update(cx_score=float(g25[~br].mean()), br_score=float(g25[br].mean()))
# v11 final: 30.0/46 and 7.0/10, the v10.1 values, so the original sentence
# ("scored similarly ... neither subset determines the pooled result") holds.
assert abs(V["br_score"] - V["cx_score"]) < 0.1, (V["br_score"], V["cx_score"])
grades = B[["mech_consistency_%s" % t for t in TAGS]].astype(str)
c4 = int((grades.iloc[:, :4].nunique(axis=1) > 1).sum())
c5 = int((grades.nunique(axis=1) > 1).sum())
stable = B.mech_consistency_threshold_stable.map(lambda x: str(x).strip().lower() == "true")
assert c5 == int((~stable).sum()), (c5, int((~stable).sum()))
V.update(c4=c4, c5=c5)
EVV = es["views"]
wv, wk, q = EVV["without_whole_variant_negatives"], EVV["without_weak_evidence_negatives"], EVV["every_expectation_E1_E2"]
assert (wv["n"], wk["n"], q["n"], q["non_disruptive"], len(wk["removed"])) == (52, 49, 24, 4, 7)
assert q["strongest_constant_labels"] == ["Both fold + PPI destabilization"] and q["margin"] < 0, q
V.update(wv=wv, wk=wk, q=q)
RSN = J("COMAVI_rubric_sensitivity.json")["reference_2_5"]
assert RSN["score"] > V["score"] and RSN["strongest_constant"] == best, RSN
src_ = ccc["structure_source"]["mechanism_pattern_score_at_2.5"]
dif = ccc["structure_source"]["difference_experimental_minus_af3"]
V.update(ex=src_["experimental"], af=src_["AlphaFold 3"], dif=dif)

# ------------------------------------------------------------------ top 20
tg = rgo["grade_breakdown_in_top_k"]
rk = pv.sort_values("isds_v1", ascending=False).reset_index(drop=True)
rk["rank"] = rk.index + 1
truth = rk.structural_ground_truth.astype(int).astype(bool)
top = rk[(rk["rank"] <= 20) & truth]
gi = {v: lab.loc[(s, v), "mech_consistency_t25"] for s, v in zip(top.system, top.variant)}
assert (sum(x == "consistent" for x in gi.values()), sum(x == "partial" for x in gi.values()),
        sum(x == "inconsistent" for x in gi.values())) == (tg["consistent"], tg["partial"], tg["inconsistent"])
inc = sorted(v for v, x in gi.items() if x == "inconsistent")
# v11 final: V11G, L749P and D351H are consistent under evidence resolution.
assert inc == sorted(["E542K", "Q542L", "N102T", "R755W", "E545K"]), inc
assert sorted(v for v, x in gi.items() if x == "consistent" and v in ("V11G", "L749P", "D351H")) == ["D351H", "L749P", "V11G"]
out4 = rk[(rk["rank"] > 20) & truth]
assert sorted(out4.variant) == sorted(["G12D", "G12V", "H718Y", "I500V"]) and \
    (int(out4["rank"].min()), int(out4["rank"].max())) == (32, 38), out4[["variant", "rank"]].values.tolist()
# The four outside the top 20: the two KRAS losses have measured magnitudes below
# the reference; H718Y and I500V have no measured energy in any committed record.
assert _mv["G12D"]["measured_kcal_mol"] < 2.5 and _mv["G12V"]["measured_kcal_mol"] < 2.5
assert _mv["H718Y"]["measured_kcal_mol"] is None and _mv["I500V"]["measured_kcal_mol"] is None


def axes_of(v):
    r = one(B, variant=v)
    return tuple(ac.MECH_LABEL_ASSERTIONS[r.mech_t25][:3]), r.mech_t25


a, l = axes_of("E542K")
assert a == (False, True, False), ("E542K", l)
assert one(B, variant="E542K").expected_mech_class == "ppi_destab_mechanism"
assert str(one(B, variant="E542K").assembled_expectation_scope) != ac.EITHER_ASSEMBLED_AXIS
for v in ("Q542L", "N102T"):
    assert axes_of(v)[1] == "Interface variant (DDG neutral)", (v, axes_of(v))
    assert one(B, variant=v).expected_mech_class in ("ppi_destab_mechanism", "mixed_structural"), v
V["top"] = tg
V["cons"] = int((g25 == 1.0).sum())
V["part"] = int((g25 == 0.5).sum())

# ------------------------------------------------------------------ text
n, k = V["n"], V["k"]
E = []
sub = lambda o, nw: E.append(("sub", o, nw))
line = lambda p, nw: E.append(("line", p, nw))

# Abstract (D1, #15, #3)
sub("We evaluated both on a literature-curated benchmark of 61 variants in 14 protein systems",
    "We evaluated both on a literature-curated benchmark of %d variants in %d protein systems" % (n, V["systems"]))
# Framing decision (11 Nov): "localizes what it detects; misses are the limit".
sub("The mechanism-pattern score was 0.661 against 0.482 for the strongest trivial strategy, but the axis the literature "
    "implicates reached the threshold for only 15 of the 29: COMAVI detects structural disruption more reliably than it "
    "localizes it. Both held across the prespecified threshold range.",
    "The mechanism-pattern score was %s against %s for the best structure-blind strategy, and exceeded it at every "
    "prespecified threshold. For %d of the %d detected mechanisms the expected axis reached the threshold; the main "
    "limitation is %d missed mechanisms, %d with measured effects below the 2.5 kcal/mol threshold."
    % (f3(V["score"]), f3(best), V["att"], V["det"], MS["n_missed"], MS["n_measured_below_reference"]))
sub("uses FoldX to score each amino acid substitution twice", "uses FoldX to score each substitution twice")
sub("Use the priority score to choose variants, then the mechanism profile to choose a stability, assembly or binding experiment.",
    "The priority score ranks variants for follow-up, and the mechanism profile indicates a stability, assembly or binding experiment.")

# Author summary
# Author summary. PLOS Computational Biology asks for 150-200 words; the v10.1
# summary was 230, so it is rewritten whole (same content, every count computed).
AS_TEXT = (
    "Genetic testing finds far more changes to single protein building blocks than any laboratory can study. Computer "
    "tools estimate which changes cause disease but rarely say how, which is what decides the next experiment. Most "
    "proteins work while bound to partners. A genetic change can make the protein itself unstable, or leave it folded "
    "but weaken its grip on a partner, and the two failures call for different experiments. Our method scores each "
    "change twice, in the protein alone and in the assembled group, to tell them apart. We tested it on %d changes in "
    "%d protein systems whose effects published experiments had already established. Judging the protein alone found "
    "damage in %d of %d disruptive changes; judging the assembled group found %d. Its predictions fully matched "
    "the published evidence for %d of the %d changes and partly for %d more, and its ranking placed %d of %d truly "
    "disruptive changes in the top 20 of %d. When it found damage, it named the right problem in %d of %d cases; its "
    "main weakness was missing %d disruptive changes, %s with only small measured effects. We describe how to choose "
    "experiments from its output and when to distrust it."
    % (n, V["systems"], V["det_iso"], V["mech"], V["det"], V["cons"], n, V["part"],
       rgo["recovery"]["structural_in_top_k"], rgo["population"]["n_structural_truth"], len(pv),
       V["att"], V["det"], MS["n_missed"], W[MS["n_measured_below_reference"]]))
assert 150 <= len(AS_TEXT.split()) <= 200, len(AS_TEXT.split())
line("Genetic testing finds far more changes to single protein building blocks", AS_TEXT)

# Introduction (#1, D2b; #2; #4)
sub("AlphaFold 3 and the Boltz models now predict complex structures for most pairs of interest, which makes this comparison "
    "possible at scale for complexes that have not been solved experimentally [@abramson2024; @wohlwend2024; @passaro2025]; "
    "whether a predicted model gives the same answer as an experimental one is itself a result, reported below.",
    "AlphaFold 3 and the Boltz models can now predict complex structures at scale, which extends structural evaluation of "
    "variants to complexes that have not been solved experimentally [@abramson2024; @wohlwend2024; @passaro2025].")
line("Several tools place variants on experimental or predicted structures",
     "Several tools place variants on experimental or predicted structures [@arnaudi2025; @khanna2021; @woodard2021; "
     "@philipp2024; @moth2026; @kwon2024; @stephenson2024], and two also use the assembled complex. ELASPIC assigns each "
     "substitution to a domain core or to an interface and predicts the corresponding change, in folding or in binding "
     "energy [@witvliet2016]. MAVISp estimates folding stability on a structure or conformational ensemble of the target "
     "protein and, in a separate module, binding energy changes for sites within 10 \u00c5 of an interface [@arnaudi2026]. "
     "Neither reports the stability of the same "
     "substitution both in the isolated subunit and in the assembled complex. That comparison, alongside the change in "
     "binding, is what COMAVI adds, and each of its outcomes points to a different experiment. A substitution that "
     "destabilizes the isolated subunit predicts a protein that folds poorly and is depleted, which calls for a stability or "
     "turnover assay. One that leaves the subunit intact but weakens an interaction predicts a protein that is present but "
     "fails to engage its partner, which calls for an affinity measurement against that partner. One that destabilizes only "
     "the assembled complex predicts a subunit that folds on its own but forms a less stable complex, which calls for a "
     "measurement on the assembled species. A single energy cannot separate these cases, because a residue buried in an "
     "interface can destabilize the fold without weakening binding. Stability and interface predictors have been benchmarked "
     "against measured energies and deep mutational scans [@buss2018; @rodrigues2019; @montanucci2019; @livesey2020; "
     "@blaabjerg2023; @sora2023; @degn2025; @vliora2026], but such benchmarks score energies, not the identification of the "
     "disrupted process.")
sub("To test both, we curated a benchmark of 61 variants in 14 protein systems",
    "To test both, we curated a benchmark of %d variants in %d protein systems" % (n, V["systems"]))
sub("Finally, a structural score invites reading as a pathogenicity prediction: on the same variants, what do the two actually measure?",
    "Finally, do structural disruption and pathogenicity coincide on the same variants?")

# Results: benchmark, populations, evidence (#3, #6-#13, D3)
sub("The benchmark contains 61 missense variants in 14 protein systems (Table 1): 49 variants in 13 multi-chain complexes "
    "and 12 in the tandem BRCA1 C-terminal (BRCT) domain",
    "The benchmark contains %d missense variants in %d protein systems (Table 1): %d variants in %d multi-chain complexes "
    "and %d in the tandem BRCA1 C-terminal (BRCT) domain" % (n, V["systems"], V["n_cx"], V["cx_systems"], V["n_br"]))
sub("Each of the 99 committed expectations is recorded with its evidence class",
    "Each of the %d committed expectations is recorded with its evidence class" % V["ledger_n"])
line("@table table1 |",
     "@table table1 | Table 1 | Benchmark composition. | Every variant carries a committed expectation on at least one "
     "modeled axis or, for four non-disruptive variants whose mechanism lies outside those axes, only at the whole-variant "
     "level, where the expectation is that no modeled axis is disrupted. Mechanisms counts the variants with a committed "
     "structural mechanism; the others are non-disruptive. Pathogenic counts the variants classified as pathogenic or likely "
     "pathogenic, by a ClinVar expert panel where one has reviewed the variant [@landrum2018] and otherwise by Franklin "
     "[@franklin]. That classification is independent of the mechanism ground truth and "
     "enters no mechanism metric; the BRCT variants carry none (\u2013). Systems are listed alphabetically. BRCT, BRCA1 "
     "C-terminal domain; VHL, von Hippel\u2013Lindau tumour suppressor; VWF, von Willebrand factor; GPIb\u03b1, platelet "
     "glycoprotein Ib \u03b1 chain. Experimental structures are identified by their Protein Data Bank (PDB) entries "
     "[@williams2001; @brzovic2001; @wu2009; @fermi1984; @tran2021; @dumas2004]; the other systems used AlphaFold 3 models.")
line("Different claims rest on different subsets of the cohort",
     "Every benchmark variant is graded. Of the %d, %d carry a committed structural mechanism and %d are curated to carry no "
     "disruption on the modeled axes, which we call non-disruptive variants. The %d variants in complexes (%d structural "
     "mechanisms and %d non-disruptive variants) also carry a structural-context tier and form the prioritization "
     "population. The %d BRCT variants (%d fold mechanisms and %d non-disruptive variants) carry none, because a single chain "
     "has no partner. The pathogenicity comparison uses the %d variants of the prioritization population with a pathogenic "
     "or benign clinical classification; the other %d carry an uncertain-significance classification or none (Fig 2A). "
     "Counts below refer to these populations by name."
     % (n, V["mech"], V["nondis"], V["n_cx"], V["cx_mech"], V["cx_nondis"], V["n_br"], V["br_mech"], V["br_nondis"],
        V["n_am"], V["not_compared"]))
line("@figure Fig 2 |",
     "@figure Fig 2 | The benchmark, and what scoring the assembled complex adds. | (A) The %d benchmark variants and the "
     "populations each analysis uses; red marks variants carrying a committed structural mechanism. (B) Detection at the "
     "2.5 kcal/mol reference by curated class, colored by the first axis to reach the threshold, taken in the order isolated "
     "subunit, assembled-complex stability, binding; for the %d non-disruptive variants the colored segments are false "
     "calls. Restricted to the %d structural mechanisms in complexes, the isolated-subunit axis detects %d and the three "
     "axes together %d. (C) Detections (above the line, of %d) and false calls (below, of %d) at each threshold, colored as "
     "in B, with false calls in lighter shades. (D) Detection, placement on the expected axis and false calls at 2.5 "
     "kcal/mol with axes removed from the analysis: the isolated-subunit axis alone, with binding added, with "
     "assembled-complex stability added, and all three. The recount uses COMAVI's own energies; it compares designs, not "
     "published tools, which use other structures and energy functions."
     % (n, V["nondis"], V["cx_mech"], V["cx_iso"], V["cx_all"], V["mech"], V["nondis"]))
E.append(("delete", "Four graded non-disruptive variants carry no committed expectation on any individual axis"))
line("The counts so far are variants.",
     "The %d variants carry %d axis expectations, one for each axis on which the literature commits a direction (S1 Table). "
     "Forty rest on quantitative measurements (%d energetic, class E1, and %d functional, E2; Methods), %d on qualitative "
     "experiments (E3), %d on population frequency (E4) and %d on inference without an axis-specific assay (E5). Evidence "
     "quality differs by direction, because a disruptive effect is far more often followed up with a dedicated biophysical "
     "measurement than a normal one is. Thus %d of the %d destabilizing expectations are quantitative, whereas all %d "
     "population and inferred expectations are neutral and %d of them belong to non-disruptive variants. We accepted these "
     "weaker classes for neutral expectations because direct measurements showing that a substitution leaves stability and "
     "binding unchanged are scarce; without them the benchmark would contain few non-disruptive variants. A non-disruptive "
     "commitment therefore records that no axis-specific assay in the primary literature supports a disruption on that axis. "
     "It does not record that a structural effect was measured and excluded. Four non-disruptive variants carry no axis "
     "expectation at all: CALM1 D96V, N98S and F142L, which act through altered calcium binding [@crotti2013; "
     "@nyegaard2012], and SMAD4 I500T, which acts through impaired monoubiquitination [@legoff2012]. They are graded on the "
     "whole-variant expectation that no modeled axis is disrupted. How far the results depend on the weaker evidence "
     "classes is reported with the mechanism-pattern score below."
     % (n, V["ledger_n"], V["e1"], V["e2"], V["e3"], V["e4"], V["e5"], V["destab_quant"], V["destab"], V["e45"], V["e45_nd"]))

# Results: what the assembled complex adds (#14; axis-design comparison)
line("COMAVI's premise is that one calculation on the isolated subunit is not enough",
     "COMAVI's premise is that one calculation on the isolated subunit is not enough, so we first measured what the assembled "
     "complex adds on the same variants at the same threshold (Fig 2B). At the 2.5 kcal/mol reference, the isolated-subunit "
     "axis alone reached the threshold for %d of the %d structural mechanisms and the three axes together for %d, at a cost "
     "of %s further false calls among the %d non-disruptive variants. Resampling whole systems puts a 95%% interval of %s to "
     "%s on the gain in detection rate, and of %s to %s once the added false calls are subtracted. The added detections came "
     "mainly from interface and mixed mechanisms, and from %d of the %d complexes rather than one or two. At every threshold "
     "the assembled axes find mechanisms that the isolated subunit misses, from %d of %d detections at 1.0 kcal/mol to %d of "
     "%d at the prediction-interval limit, and they account for every additional false call, a cost that falls from %d "
     "variants at 1.0 kcal/mol to %d at the limit (Fig 2C; S3 Table, threshold_sweep sheet)."
     % (V["det_iso"], V["mech"], V["det"], W[V["added_fc"]], V["nondis"],
        f2(SC["detection_gain"]["ci95"][0]), f2(SC["detection_gain"]["ci95"][1]),
        f2(SC["net_gain"]["ci95"][0]), f2(SC["net_gain"]["ci95"][1]), V["n_sys_added"], V["cx_systems"],
        asm["t10"], sw["t10"]["detection_k"], asm["tSAP"], sw["tSAP"]["detection_k"], V["fc10"], V["fcS"])
     + "\n\n%%from: new (v11; axis-design comparison, author decision)\n"
     "Of the two assembled axes, assembled-complex stability is the one that neither ELASPIC nor MAVISp computes. Removing it "
     "from the analysis leaves isolated-subunit stability and binding, the two quantities MAVISp's stability and "
     "local-interaction modules estimate. Without it, detection falls from %d to %d mechanisms and the number placed on the "
     "expected axis from %d to %d, while %s false calls are avoided (Fig 2D). Of the %s mechanisms detected only through "
     "assembled-complex stability, %s are placed on the expected axis; the fifth, PIK3CA E542K, is an interface "
     "destabilizer that it reports as instability of the assembly rather than as loss of binding."
     % (V["d_all"], V["d_ib"], V["a_all"], V["a_ib"], W[V["fc_saved"]], W[n_only], W[V["only_att"]]))
assert n_only == 5

# Results: mechanism-pattern score (#15; framing decision: margin kept as headline)
sub("## COMAVI recovers literature-defined mechanisms but localizes them less reliably than it detects them",
    "## COMAVI localizes most of the mechanisms it detects")
line("Each graded variant received a three-level grade",
     "Each graded variant received a three-level grade (Table 3): 1 when COMAVI recovered the expected combination of "
     "effects; 0.5 when it recovered part of it, missing one of several expected effects or adding one the literature does "
     "not support; and 0 when it recovered none of them or predicted a change in the opposite direction. Table 2 defines "
     "this mechanism-pattern score and every other metric reported below. At the 2.5 kcal/mol (10.5 kJ/mol) reference "
     "threshold, the mechanism-pattern score was %.1f of %d (%s). To judge this value, we graded structure-blind strategies, "
     "which ignore the structure and give every variant the same mechanism label, for each of the %d labels COMAVI can "
     "assign. The strongest, calling no variant disrupted, scores %s by accepting all non-disruptive variants and missing "
     "every mechanism; the three labels that name a structural or interface site without an energy change score the same. "
     "Calling every variant destabilized on all three axes scores %s, and permuting the ground truth gives a mean of %s "
     "(standard deviation %s; p = %.4f; S2 Text). COMAVI exceeds the best structure-blind strategy by %s. Because that "
     "strategy scores exactly the non-disruptive fraction of the cohort, this margin, rather than the raw score, is the "
     "quantity that can be compared across cohorts of different composition. The margin is positive at every threshold and "
     "with any one system left out (%s to %s), but resampling whole systems puts a 95%% interval of %s to %s on it, and %d%% "
     "of resamples give a margin of zero or less. Of the %.1f grade points by which COMAVI exceeds that strategy, the BRCT "
     "domain contributes %.1f and the hemoglobin tetramer %.1f, while %s systems reduce it (S1 Text). The benchmark "
     "therefore shows the advantage on these 14 systems but cannot establish that it would hold in others."
     % (k, n, f3(V["score"]), V["n_labels"], f3(best), f3(V["all3"]), f3(V["perm_mean"]), f3(V["perm_sd"]), V["perm_p"],
        f3(V["margin"]), f3(V["lm"][0]), f3(V["lm"][1]), f3(V["m_ci"][0]), f3(V["m_ci"][1]), round(100 * V["m_p"]),
        V["c_total"], V["c_brct"], V["c_hb"], W[V["n_neg"]]))
line("The pooled score combines two capabilities that behave differently.",
     "The pooled score combines two capabilities that behave differently. The grade averaged %s over the %d "
     "structural-mechanism variants and %s over the %d non-disruptive variants. Asked as binary questions, the same calls "
     "give three rates. At least one axis reached the threshold for %d of the %d mechanisms (detection, %s), and the axis "
     "the literature implicates reached it for %d (attribution, %s). No axis reached it for %d of the %d non-disruptive "
     "variants (correct rejection, %s). The graded structural arm lies below binary attribution at every uniform threshold, "
     "by %s to %s, so partial credit does not inflate it. COMAVI therefore localizes most of what it detects, and its "
     "principal limitation is the %d mechanisms that no axis detects. For %s of these the measured effect is itself below "
     "the 2.5 kcal/mol threshold (%s to %s kcal/mol: KRAS G12D and G12V, which weaken binding to RAF, and BRCA1 R1751Q and "
     "L1664P, which destabilize the BRCT fold), so even an exact prediction would not reach it. For %s others the measured "
     "effect is large but the predicted energy falls short on every axis (hemoglobin N102T, %s kcal/mol measured, and BRCA1 "
     "V1736A, %s kcal/mol), and %s have no measured energy (S3 Table, missed_mechanisms sheet)."
     % (f3(arms["t25"][0]), V["mech"], f3(arms["t25"][1]), V["nondis"], V["det"], V["mech"], f3(V["det"] / V["mech"]),
        V["att"], f3(V["att"] / V["mech"]), V["rej"], V["nondis"], f3(V["rej"] / V["nondis"]),
        f3(min(gaps.values())), f3(max(gaps.values())), MS["n_missed"], W[len(_below)],
        f2(V["below_rng"][0]), f2(V["below_rng"][1]), W[len(_large)],
        ("%.1f" % _mv["N102T"]["measured_kcal_mol"]), ("%.1f" % _mv["V1736A"]["measured_kcal_mol"]), W[len(_none)]))
line("The localization errors concentrate in one class.",
     "Localization errors are few. Of the %d detected mechanisms, %s reached the threshold only on an axis other than the "
     "one the literature implicates: PIK3CA E542K, an interface destabilizer reported as instability of the assembly, and "
     "MLH1 R755W, a fold mechanism reported as strengthened binding (discussed with the ranking below). Resampling whole "
     "systems puts a 95%% interval of %s to %s on this wrong-axis fraction of the %d mechanisms. Among the %d interface "
     "variants curated as destabilizers, the binding axis reached the threshold for %d and the assembled-complex stability "
     "axis for %d. An interface substitution can genuinely raise the folding energy of the assembled coordinates, so the two "
     "assembled axes are not independent, and for %s interface expectations the evidence cannot separate them either: they "
     "rest only on assays that report how much complex forms, such as co-immunoprecipitation, pull-down or two-hybrid "
     "assays, and a change in the amount of complex can come from a change in binding or in the stability of the "
     "assembly. These expectations "
     "are therefore met by a change in the measured direction on either assembled axis (Methods). %s of the %d localized "
     "mechanisms, BRCA1 V11G, MLH1 L749P and SMAD4 D351H, are credited in this way; without the rule they would count as "
     "detected on the wrong axis (S2 Text)."
     % (V["det"], W[len(V["wrong_axis"])], f2(V["loc_gap_ci"][0]), f2(V["loc_gap_ci"][1]), V["mech"],
        V["n_if_destab"], int(ifd.bind_fires.sum()), int(ifd.fold_fires.sum()), W[V["n_scope"]],
        Wc[len(V["credited_by_rule"])], V["att"]))
# The subset sentence stays as written: v11 final reproduces its values.
assert (V["n_cx"], f3(V["cx_score"]), V["n_br"], f3(V["br_score"])) == (46, "0.652", 10, "0.700")
_ax = lambda a: (t25["axis_%s_k" % a], t25["axis_%s_n" % a])
_rates = [t25["axis_%s_k" % a] / t25["axis_%s_n" % a] for a in ("monomer", "fold", "binding")]
_energ = [sum(t25["axis_%s_k" % a] for a in ac.ENERGETIC_ROWS), sum(t25["axis_%s_n" % a] for a in ac.ENERGETIC_ROWS)]
assert set(ac.ENERGETIC_ROWS) == {"monomer", "fold", "binding", "assembled_either"}
sub("Direction agreement was similar on the three axes: 19 of 26 expectations for isolated-subunit stability, 17 of 22 for "
    "assembled-complex stability and 20 of 29 for binding (56 of the 77 axes overall that pass confidence gating in a "
    "graded variant, counting axes rather than variants), rates spanning 0.083.",
    "Direction agreement was similar on the three axes: %d of %d expectations for isolated-subunit stability, %d of %d for "
    "assembled-complex stability and %d of %d for binding, rates spanning %s. The %d expectations measured only as an "
    "amount of complex were scored jointly, and a call in the measured direction on either assembled axis met %d of them "
    "(%d of the %d axes overall that pass confidence gating in a graded variant, counting axes rather than variants)."
    % (_ax("monomer") + _ax("fold") + _ax("binding") + (f3(max(_rates) - min(_rates)), _ax("assembled_either")[1],
       _ax("assembled_either")[0], _energ[0], _energ[1])))
wv_k = round(wv["comavi_score"] * wv["n"] * 2) / 2
line("No single protein system carries the result.",
     "No single protein system carries the result. Resampling whole systems gave a 95%% interval of %s\u2013%s for the "
     "mechanism-pattern score, and leaving out one system at a time gave values from %s to %s (Fig 3). In a sensitivity "
     "analysis that gives half credit to calls whose asserted changes all fall on axes the literature never tested, the "
     "score rises to %s (S2 Text). Without the four "
     "whole-variant non-disruptive variants, the score is %.1f of %d (%s), against %s for the strongest structure-blind "
     "strategy on the same variants. Without all %s non-disruptive variants that rest only on population, inferred or "
     "whole-variant evidence, the margin narrows to %s (95%% interval %s to %s). On the %d variants whose expectations are "
     "all quantitative, only %s of them non-disruptive, calling every variant destabilized on all three axes scores %s, "
     "above COMAVI's %s. Because these subsets differ widely in class balance, S1 Text also reports the balanced score, the "
     "mean of the two arms, on which no structure-blind strategy can exceed 0.5: COMAVI scores %s without the weakly "
     "evidenced non-disruptive variants and %s on the quantitative subset. Because the %d variants come from 14 clustered "
     "systems, %s measures agreement on this benchmark, not general accuracy."
     % (f3(V["ci"][0]), f3(V["ci"][1]), f3(V["loso"][0]), f3(V["loso"][1]), f3(RSN["score"]), wv_k, wv["n"],
        f3(wv["comavi_score"]),
        f3(wv["strongest_constant"]), W[len(wk["removed"])], f3(wk["margin"]), f3(wk["margin_ci95"][0]),
        f3(wk["margin_ci95"][1]), q["n"], W[q["non_disruptive"]], f3(q["strongest_constant"]), f3(q["comavi_score"]),
        f3(wk["comavi_balanced"]), f3(q["comavi_balanced"]), n, f3(V["score"])))
sub("the vertical line marks the observed 0.661", "the vertical line marks the observed %s" % f3(V["score"]))
sub("dashed lines mark the 95% interval (0.573\u20130.755)",
    "dashed lines mark the 95%% interval (%s\u2013%s)" % (f3(V["ci"][0]), f3(V["ci"][1])))
sub("the mechanism-pattern score was 20.0 of 28 (0.714) on the experimental systems and 17.0 of 28 (0.607) on the predicted ones",
    "the mechanism-pattern score was %.1f of %d (%s) on the experimental systems and %.1f of %d (%s) on the predicted ones"
    % (V["ex"]["score_total"], V["ex"]["variants"], f3(V["ex"]["score"]), V["af"]["score_total"], V["af"]["variants"],
       f3(V["af"]["score"])))
sub("puts a 95% interval of \u22120.080 to 0.257 on the difference of 0.107",
    "puts a 95%% interval of %s to %s on the difference of %s"
    % (f3(dif["ci95_system_resampling"][0]), f3(dif["ci95_system_resampling"][1]), f3(dif["value"])))

# Results: stringency
sub("the mechanism-pattern score was 0.562, 0.616, 0.652, 0.661 and 0.652.",
    "the mechanism-pattern score was %s." % (", ".join(f3(sw[t]["whole_variant"]) for t in TAGS[:4]) + " and "
                                             + f3(sw["tSAP"]["whole_variant"])))
_cb = J("COMAVI_threshold_sweep.json")["cluster_bands"]
assert all(max(_cb[a]["lo"], _cb[b]["lo"]) <= min(_cb[a]["hi"], _cb[b]["hi"])
           for i, a in enumerate(TAGS) for b in TAGS[i + 1:]), "system-resampling intervals overlap at every threshold"
sub("## Stringency costs attribution, not detection",
    "## Stringency trades detection for correct rejection")
_seq = lambda xs: ", ".join(str(x) for x in xs[1:-1]) + " and " + str(xs[-1])
line("Raising the threshold changed the three binary rates in different ways.",
     "Raising the threshold changed the three binary rates in different ways. Detection fell: at least one axis reached the "
     "threshold for %d of the %d mechanisms at 1.0 kcal/mol, for %d at 1.5, %d at 2.0, %d at 2.5 and %d at the "
     "prediction-interval limit. Attribution fell less, from %d to %s, so the number of mechanisms detected only on another "
     "axis shrank from %d to %d. Correct rejection rose at every step, from %d of %d to %s. Raising the threshold therefore "
     "gains correct rejections by missing disruptions, not by placing them on the wrong axis (Fig 4B). The structural arm "
     "fell from %s to %s and the non-disruptive arm rose from %s to %s, and the two cross between %s (Fig 4A)."
     % (V["det_t"][0], V["mech"], V["det_t"][1], V["det_t"][2], V["det_t"][3], V["det_t"][4], V["att_t"][0],
        _seq(V["att_t"]), V["wa_t"][0], V["wa_t"][-1], V["rej_t"][0], V["nondis"], _seq(V["rej_t"]),
        f3(arms["t10"][0]), f3(arms["tSAP"][0]), f3(arms["t10"][1]), f3(arms["tSAP"][1]), V["cross"]))
_bk = [sw[t]["axis_binding_k"] for t in TAGS]
_fk = [sw[t]["axis_fold_k"] for t in TAGS]
_mk = [sw[t]["axis_monomer_k"] for t in TAGS]
assert len(set(sw[t]["axis_binding_n"] for t in TAGS)) == 1 and len(set(_bk)) == 1, _bk
assert max(_fk) == sw["t25"]["axis_fold_k"] and _fk[0] == min(_fk), _fk
_mbest = [TLAB[t] for t in TAGS if sw[t]["axis_monomer_k"] == max(_mk)]
assert _mbest == ["2.0", "2.5"], _mbest
line("The three axes prefer different thresholds (Fig 4C).",
     "The two stability axes prefer stringent thresholds, and binding is indifferent (Fig 4C). Assembled-complex agreement "
     "rose from %d of %d at 1.0 kcal/mol to %d at the reference, and isolated-subunit agreement peaked at 2.0 and 2.5 "
     "kcal/mol (%d of %d). Binding agreement was %d of %d at every threshold. The denominators do not change across "
     "thresholds. Per-axis thresholds therefore offer little: holding the two stability axes at the reference and lowering "
     "the binding threshold anywhere in the prespecified range leaves attribution at %d of %d and interface attribution at "
     "%d of %d, and at 1.0 or 1.5 kcal/mol costs %s correct rejections. Each of the %d interface mechanisms not localized at "
     "the reference carries a predicted binding energy below 1.0 kcal/mol, the most permissive threshold in the prespecified "
     "range (median %s), so no placement of a binding threshold reaches them. For those variants the predicted binding "
     "energy is small where the literature reports a changed interaction, which locates the shortfall in the energy "
     "prediction rather than in the calling layer."
     % (_fk[0], sw["t25"]["axis_fold_n"], _fk[3], max(_mk), sw["t25"]["axis_monomer_n"], _bk[0],
        sw["t25"]["axis_binding_n"], V["att"], V["mech"], V["pa_if"][0], V["pa_if"][1], W[V["pa_lost"]], V["pa_ni"],
        f2(V["pa_med"])))

# Results: top 20
line("A high rank does not guarantee a correct mechanism.",
     "A high rank does not guarantee a correct mechanism. Of the 17 structural mechanisms in the top 20, %d are graded "
     "consistent, %d partial and %d inconsistent. The four outside it, at ranks 32 to 38, are KRAS G12D and G12V, whose "
     "measured losses of RAF binding (%s and %s kcal/mol) are below the 2.5 kcal/mol reference, and MLH1 H718Y and SMAD4 "
     "I500V, whose effects were measured without an energy. Three of the five inconsistent grades are expected losses of "
     "binding that COMAVI either placed on another axis (PIK3CA E542K, as instability of the assembly) or called at an "
     "interface site with no energy change (MLH1\u2013PMS2 Q542L and hemoglobin N102T). "
     "The other two, MLH1 R755W and PIK3CA E545K, fail in a different way. FoldX returns a large favorable interaction "
     "energy for each (\u22129.7 and \u22124.6 kcal/mol), which COMAVI reads as strengthened binding. For E545K the axis is "
     "right but the direction is wrong, because releasing the inhibitory contact with p85 is energetically favorable "
     "[@zhao2008]. For R755W the literature reports reduced abundance with binding preserved [@kosinski2010; "
     "@larsenledet2025], and neither stability prediction reaches the threshold (0.97 kcal/mol for the isolated subunit and "
     "\u22120.94 in the complex). These two variants would send a user to the right protein with the wrong experiment "
     "(S3 Text)." % (V["top"]["consistent"], V["top"]["partial"], V["top"]["inconsistent"],
                     "%.2f" % _mv["G12D"]["measured_kcal_mol"], "%.2f" % _mv["G12V"]["measured_kcal_mol"]))
# The FoldX energies quoted for R755W and E545K, from the canonical.
_fc_cols, _bc_cols = bda.axis_columns(can, "fold"), bda.axis_columns(can, "binding")
_e = lambda s, v: (lambda c: (round(float(c.ddg_monomer), 2), round(float(bda.strongest(c, _fc_cols, "ddg_fold_")), 2),
                              round(float(bda.strongest(c, _bc_cols, "ddg_binding_")), 1)))(
    can.set_index(["system", "variant"]).loc[(s, v)])
assert _e("mlh1_pms2", "R755W") == (0.97, -0.94, -9.7) and _e("pi3k", "E545K")[2] == -4.6, (_e("mlh1_pms2", "R755W"),
                                                                                         _e("pi3k", "E545K"))

# Results: pathogenicity (v11 rule: ClinVar expert panel, then Franklin)
line("Clinical annotation and modeled structural disruption come apart inside this cohort.",
     "Clinical annotation and modeled structural disruption come apart inside this cohort. Of the %d variants classified as "
     "pathogenic among the %d in the comparison, %d are curated to carry no disruption on any modeled axis "
     "and %d carry a structural mechanism. Half of the pathogenic class is therefore invisible to a structural method by "
     "construction, which is the strongest statement the benchmark makes about the two properties being distinct."
     % (V["P"], V["n_am"], V["P_nondis"], V["P_struct"]))
line("Rank-based comparison with a pathogenicity predictor is correspondingly limited.",
     "Rank-based comparison with a pathogenicity predictor is correspondingly limited. %s of the %d benign variants in this "
     "set entered the benchmark as curated non-disruptive controls, so any structural score separates them from pathogenic "
     "structural mechanisms partly by construction, and a pathogenicity AUC computed here measures the cohort as much as "
     "the method. The exception is informative: MLH1 H718Y is clinically benign and carries a committed monomer-fold "
     "destabilization, because the source reports reduced thermal stability and a half-life of 64%% of wild type alongside "
     "proficient mismatch repair. A variant can be structurally destabilized and clinically neutral at once, which is the "
     "separation the two outputs are meant to preserve. With that caveat, AlphaMissense ranked pathogenicity better (ROC "
     "AUC %s) than the energy component (%s), the priority score (%s) or the context component (%s) (Fig 7). This set is "
     "the part of the prioritization population that carries a pathogenic or benign classification, but the positive "
     "classes differ, clinical label against curated mechanism, so these values are not comparable with the structural AUC "
     "of %s." % (Wc[V["B_nondis"]], V["Bn"], f3(AUC["am"]), f3(AUC["en"]), f3(AUC["pr"]), f3(AUC["cx"]),
                  f3(pm_isds["roc_auc"])))
sub("for the 47 complex variants with a clinical label and an AlphaMissense score",
    "for the %d variants in complexes with a pathogenic or benign clinical classification" % V["n_am"])
av = V["am_vals"]
line("Discordant variants show what each view adds.",
     "Discordant variants show what each view adds. HBB E6V is pathogenic but scores low on both AlphaMissense (%s) and the "
     "priority score (%s), because its mechanism is absent from the supplied model. CFH R78G has a high priority score (%s) "
     "and a named binding disruption. PIK3CA H1047R sits in Tier 2 context but has no energy call on any axis at the "
     "reference threshold, so COMAVI ranks it on context alone and names no mechanism. Both of these clinically pathogenic "
     "variants, R78G and H1047R, fall in the AlphaMissense ambiguous range (%s and %s), as does one benign variant, MSH2 "
     "N127S. Where the pathogenicity score is least decided, the structural profile supplies a named mechanism for one "
     "variant and a place to look for another."
     % (f3(av["E6V"][0]), f3(av["E6V"][1]), f3(av["R78G"][1]), f3(av["R78G"][0]), f3(av["H1047R"][0])))

# Discussion and conclusion
sub("17 of the 56 graded variants change grade across the four uniform thresholds, and 20 across all five operating points",
    "%d of the %d variants change grade across the four uniform thresholds, and %d across all five operating points"
    % (V["c4"], n, V["c5"]))
assert V["det"] / V["mech"] > 0.5
sub("On this benchmark it recovered most literature-defined mechanisms, ranked structural mechanisms near the top, and "
    "failed most often by placing a detected disruption on the wrong axis.",
    "On this benchmark it detected most literature-defined mechanisms and localized most of those it detected, ranked "
    "structural mechanisms near the top, and failed most often by missing a disruption whose predicted energy stayed "
    "below the threshold.")
assert MS["n_missed"] > len(V["wrong_axis"]) and MS["n_missed"] > V["nondis"] - V["rej"], "missing is the commonest failure"
# Discussion: MAVISp (checked against the published paper, Protein Science 2026, PMC13055194).
sub("whereas the MAVISp stability module uses one structure per target and its interaction module treats interfaces "
    "separately.",
    "whereas MAVISp estimates folding stability on a structure or conformational ensemble of the target protein and, in "
    "its local-interaction module, binding free energy for sites within 10 \u00c5 of an interface.")
sub("MAVISp is validated against functional assays and clinical labels, and its authors note that interaction data are too "
    "scarce and heterogeneous to benchmark its binding module.",
    "MAVISp is validated against functional assays and clinical labels; its authors set binding-energy thresholds from "
    "the expected error of the predictors and note that measured effects on protein\u2013protein interactions are scarce "
    "and heterogeneous enough to complicate reliable benchmarking.")
sub("whereas MAVISp is a curated database to which new proteins are added.",
    "whereas MAVISp is a curated database to which trained biocurators add proteins, including targets that users request.")
# Limitations: sub-threshold expectations, at full precision from the generator
# (expectations measured only as an amount of complex are scored jointly there too).
_sb_spec = importlib.util.spec_from_file_location("build_subthreshold_evidence",
                                                  REPO / "scripts" / "build_subthreshold_evidence.py")
_sb_mod = importlib.util.module_from_spec(_sb_spec)
_sb_spec.loader.exec_module(_sb_mod)
SB = _sb_mod.build(digits=None)["axis_level"]
_sb_rec = J("COMAVI_subthreshold_evidence.json")["axis_level"]
assert all(round(SB[k], 3) == _sb_rec[k] for k in ("auc_within_band", "p_auc_at_or_below_chance")) and \
    [round(x, 3) for x in SB["auc_within_band_ci95"]] == _sb_rec["auc_within_band_ci95"], "full precision vs record"
assert SB["auc_within_band_ci95"][0] <= 0.5 <= SB["auc_within_band_ci95"][1], "interval includes chance"
sub("Twelve committed destabilizing axes fall below even 1.0 kcal/mol, and within that range the predicted magnitude "
    "does not distinguish them from axes the literature calls intact (AUC 0.63, 95% confidence interval 0.39 to 0.81).",
    "%s committed destabilizing expectations fall below even 1.0 kcal/mol, and within that range the predicted magnitude "
    "does not distinguish them from those the literature calls intact (AUC %s, 95%% confidence interval %s to %s)."
    % (Wc[SB["below_and_committed_destabilizing"]], f2(SB["auc_within_band"]), f2(SB["auc_within_band_ci95"][0]),
       f2(SB["auc_within_band_ci95"][1])))
# Fig 4 caption: panel C plots the three axis-specific rates.
sub("(C) Direction agreement on each axis (26 isolated-subunit, 22 assembled-complex and 29 binding expectations); open "
    "circles mark the best value for each axis.",
    "(C) Direction agreement on each axis (%d isolated-subunit, %d assembled-complex and %d binding expectations; the %d "
    "expectations measured only as an amount of complex are scored jointly and given in S3 Table); open circles mark the "
    "best value for each axis." % (sw["t25"]["axis_monomer_n"], sw["t25"]["axis_fold_n"], sw["t25"]["axis_binding_n"],
                                   sw["t25"]["axis_assembled_either_n"]))
sub("the interface-class localization error, interface ablation of the tier",
    "localization by mechanism class and under the evidence-resolution rule, the untested-axis sensitivity analysis, "
    "interface ablation of the tier")

# ------------------------------------------------------------------ clinical-label rule (Methods)
CV = pd.read_csv(REPO / "data" / "analysis_inputs" / "clinvar" / "COMAVI_clinvar_records.csv", keep_default_na=False, dtype=str)
assert set(zip(CV.system, CV.variant)) == set(zip(B[~br].system, B[~br].variant)), "records cover every complex variant"
_dates = sorted(CV.retrieved.unique())
assert len(_dates) == 1, _dates
V["cv_date"] = "%d %s %s" % (int(_dates[0][8:]), ["January", "February", "March", "April", "May", "June", "July", "August",
                                                    "September", "October", "November", "December"][int(_dates[0][5:7]) - 1],
                             _dates[0][:4])
_ep = CV[CV.review_status.isin(["reviewed by expert panel", "practice guideline"])]
_cls = lambda d: {"pathogenic": "pathogenic", "likely pathogenic": "pathogenic", "pathogenic/likely pathogenic": "pathogenic",
                  "benign": "benign", "likely benign": "benign", "benign/likely benign": "benign",
                  "uncertain significance": "vus"}[d.strip().lower()]
_epc = _ep.assign(c=_ep.classification.map(_cls)).groupby(["system", "variant"]).c.apply(set)
_fr = lambda x: ("pathogenic" if "pathogenic" in str(x).lower() and "benign" not in str(x).lower()
                 else "benign" if "benign" in str(x).lower() and "pathogenic" not in str(x).lower() else "other")
_frk = can.set_index(["system", "variant"]).franklin
_agree = [k for k, c in _epc.items() if len(c) == 1 and _fr(_frk.loc[k]) == next(iter(c))]
_change = [k[1] for k, c in _epc.items() if len(c) == 1 and _fr(_frk.loc[k]) != next(iter(c))]
_conflict = [k[1] for k, c in _epc.items() if len(c) > 1]
V.update(ep_n=len(_epc), ep_agree=len(_agree))
assert (len(_epc), len(_agree), _change, _conflict) == (14, 12, ["G857A"], ["M326I"]), (len(_epc), len(_agree), _change, _conflict)
assert one(B, variant="G857A").phenotype == "benign"
# The six complex-amount expectations and their assays (Methods).
_ca = led[(led.complex_amount_assay.fillna("") != "") & led.expected_token.isin(["destab", "stab"])]
assert sorted(zip(_ca.system, _ca.variant)) == sorted(zip(_scope.system, _scope.variant)), "scope = directional complex-amount rows"
assert all(re.search(r"co-IP|pull-down|two-hybrid|Y2H|TR-FRET", a) for a in _ca.complex_amount_assay), list(_ca.complex_amount_assay)

# Methods
sub("S1 Table gives the class, basis and citation of all 99 committed expectations",
    "S1 Table gives the class, basis and citation of all %d committed expectations" % V["ledger_n"])
E.append(("insert_after", "The primary literature defined the expected direction on each applicable axis",
          "%%from: new (v11; reviewer comments 3 and 6-10)\n"
          "Five scored candidates were excluded from the benchmark because no gradeable expectation could be defined for "
          "them. For TNNI3 R162W, VHL W117R and MSH2 G674R, the primary literature commits no direction on any modeled axis. "
          "BRCA1 R1699L and R1699Q disrupt phosphopeptide recognition, a protein interaction whose partner is absent from "
          "the supplied structure. Variants whose mechanism lies outside the processes COMAVI scores, such as calcium "
          "binding or monoubiquitination, were instead retained and graded on the expectation that no modeled axis is "
          "disrupted. The five are listed with their COMAVI outputs in S2 Table (excluded_variants sheet).\n\n"
          "%%from: new (v11; clinical labels and AlphaMissense provenance)\n"
          "Clinical labels enter only the pathogenicity comparison, and each follows one rule. A ClinVar classification "
          "reviewed by an expert panel or given in a practice guideline takes precedence [@landrum2018]; otherwise the "
          "Franklin classification is used [@franklin]. ClinVar records for the %d variants in complexes were retrieved "
          "through NCBI E-utilities on %s and are distributed with the code. %s variants carry an expert-panel "
          "classification: it agrees with Franklin for %s, changes PMS2 G857A from a Franklin variant of uncertain "
          "significance to benign, and for PIK3R1 M326I is split between two records for different nucleotide changes, "
          "so Franklin is used. Pathogenic or likely pathogenic variants are pathogenic, benign "
          "or likely benign variants are benign, and variants of uncertain significance or without a classification are "
          "not compared. AlphaMissense scores were "
          "taken from the published release (doi:10.5281/zenodo.8208688) [@cheng2023] at UniProt canonical positions, with "
          "the HBB variants converted from mature-chain numbering and the two reverse-named controls described below looked "
          "up as their clinical substitutions." % (V["n_cx"], V["cv_date"], Wc[V["ep_n"]], W[V["ep_agree"]])))
sub("Detection, attribution and correct rejection ask the same question as binary rates (Table 2).",
    "Detection, attribution and correct rejection ask the same question as binary rates (Table 2). Six interface "
    "expectations rest only on assays that report how much complex forms (co-immunoprecipitation, pull-down, two-hybrid "
    "or TR-FRET), and the evidence ledger names the assay for each (S1 Table). A change in the amount of complex can come "
    "from weaker binding or from a less stable assembled complex, which COMAVI computes separately, so for these six a "
    "change in the measured direction on either assembled axis meets the expectation. This holds wherever an "
    "expectation is scored: in grading, in attribution, in the sub-threshold analysis and in direction agreement, where "
    "the six form a joint row and each axis's own row counts only expectations specific to that axis. The rule is fixed by the evidence before any call is considered and applies "
    "equally to every structure-blind strategy; a call in the opposite direction on an assembled axis earns nothing from "
    "it. It is not applied in reverse: a normal amount of complex supports the committed binding expectation but cannot "
    "show that the assembled complex is stable, so no assembled-complex expectation was added for non-disruptive variants. As a sensitivity analysis only, "
    "a call that recovers no expected effect but asserts changes solely on axes the literature never tested, and "
    "contradicts nothing measured, was credited 0.5 (S2 Text).")
sub("1 for a fully consistent whole-variant pattern, 0.5 for a partly correct multi-axis pattern and 0 for an inconsistent pattern",
    "1 for a fully consistent whole-variant pattern, 0.5 for a pattern that recovers at least one expected effect but "
    "misses another or adds an unsupported one, and 0 for a pattern that recovers no expected effect or predicts the "
    "opposite direction")
sub("For 53 of the 61 that is the numbering of the UniProt reference sequence",
    "For %d of the %d that is the numbering of the UniProt reference sequence" % (V["uniprot_named"], n))

# Supporting information captions
sub("All 99 committed axis expectations, and the four non-disruptive variants graded at the whole-variant level",
    "All %d committed axis expectations of the benchmark variants, and the four non-disruptive variants graded at the "
    "whole-variant level" % V["ledger_n"])
sub("evidence composition, model scope, structural-context comparators, denominator sensitivity and fold\u2013binding "
    "disagreements, one sheet each",
    "evidence composition, model scope, the five excluded candidates, structural-context comparators and fold\u2013binding "
    "disagreements, one sheet each")
sub("FoldX commands, tier weights and cut-offs, seeds and draw counts, and the evidence-restricted sensitivity analyses.",
    "FoldX commands, tier weights and cut-offs, seeds and draw counts, the evidence-restricted sensitivity analyses, each "
    "against its own strongest structure-blind strategy, and the contribution of each system to the margin.")
sub("Trivial strategies and the permutation null",
    "Structure-blind strategies for all %d labels and the permutation null" % V["n_labels"])

# ------------------------------------------------------------------ apply
text = SRC.read_text()
lines = text.split("\n")


def find(prefix):
    hits = [i for i, l in enumerate(lines) if l.startswith(prefix)]
    assert len(hits) == 1, (prefix[:60], len(hits))
    return hits[0]


for kind, a, *b in [e for e in E if e[0] != "sub"]:
    if kind == "line":
        lines[find(a)] = b[0]
    elif kind == "delete":
        i = find(a)
        assert lines[i - 1].startswith("%from:") and lines[i - 2] == "", lines[i - 1][:40]
        del lines[i - 2:i + 1]
    elif kind == "insert_after":
        i = find(a)
        lines[i + 1:i + 1] = [""] + b[0].split("\n")
text = "\n".join(lines)
for kind, old, new in [e for e in E if e[0] == "sub"]:
    assert text.count(old) == 1, (old[:70], text.count(old))
    text = text.replace(old, new)

# Tokens that were true in v10.1 and are false in v11 final, and phrases of the
# superseded localization story. Any occurrence stops the run, except the
# excluded candidates, which are named exactly once (the Methods paragraph).
RETIRED = ["15 of the 29", "for 15 (", "0.083", "of the 77", "20 of 29", "0.562", "0.616", "0.672", "0.466",
           "1.5 and 2.0 kcal/mol", "0.573", "0.755", "0.641 to 0.679", "from 0.641 to", "0.058", "0.373", "0.464", "0.903", "0.851", "0.836",
           "0.769", "0.080 to 0.257", "4 of 12", "eight inconsistent", "more reliably than it localizes",
           "less reliably than it detects", "wrong problem", "costs attribution", "interface shortfall but",
           "most promising single change", "trivial", " 61 ", "of 61", " 47 ", " 49 ", " 99 ", "Franklin classifies",
           "32 and partly for 10", "Twelve committed destabilizing", "AUC 0.63", "Detection was nearly flat", "gives up 0.069", "(22 of 29)"]
ALLOWED_ONCE = ["R162W", "W117R", "G674R"]
body = "\n".join(l for l in text.split("\n") if not l.startswith("%"))
print("edits: %d" % len(E))
_tl = text.split("\n")
_abs_line = [l for l in _tl if l.startswith("Pathogenicity predictors rank missense variants")]
assert len(_abs_line) == 1 and len(_abs_line[0].split()) <= 300, ("abstract over 300 words", len(_abs_line[0].split()))
print("abstract words: %d | author summary words: %d" % (len(_abs_line[0].split()), len(AS_TEXT.split())))
bad = []
for t in RETIRED:
    for i in [j for j in range(len(body)) if body.startswith(t, j)][:6]:
        bad.append("  residual %-14r ...%s..." % (t, body[max(0, i - 70):i + 50].replace("\n", " ")))
for t in ALLOWED_ONCE:
    if body.count(t) != 1:
        bad.append("  %r appears %d times (expected once, in Methods)" % (t, body.count(t)))
print("\n".join(bad) if bad else "residual scan: clean")
assert not bad, "stale tokens remain"
audit = REPO / "docs" / "audit" / "COMAVI_v11_manuscript_values.json"
if ARGS.dry_run:
    preview = pathlib.Path(__file__).resolve().with_name("v11_manuscript_preview.txt")
    preview.write_text(text)
    print("dry run: repository untouched; preview written to", preview.name)
else:
    SRC.write_text(text)
    audit.write_text(json.dumps(V, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
    print("wrote", SRC.relative_to(REPO), "and", audit.relative_to(REPO))
