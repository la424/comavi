#!/usr/bin/env python3
"""Build the three supporting-information workbooks and the package manifest.

WHY THIS EXISTS
---------------
Until v11 the S1-S3 Table workbooks were assembled by hand, and the v9.2 set
had drifted from the records in ways no gate could see: S1 Table listed the
whole 99-row ledger rather than the benchmark's 97 expectations, with reference
numbers from an earlier manuscript; S2 and S3 Table described a 47-variant
prioritization population and a 57-variant grading cohort; the operating-point,
robustness, evidence-state and pathogenicity sheets carried pre-v11 values; the
threshold_sweep caption stopped mid-sentence; and the README counted eleven
items and five figures for a package of ten items and four figures. This script
builds every sheet from the committed records, so the workbooks follow the
records the way the manuscript and the supporting texts do:

  records -> submission/COMAVI_v*_SI/S{1,2,3}_Table.xlsx
          -> SI_INDEX.csv, README.txt and SHA256SUMS.txt in the same directory

Statistics are rounded once, half to even, on the exact value (grades are
multiples of 0.5, so most are exact fractions). Text that no record holds --
the model-scope rules, the interpretation column of fold_binding_disagreements,
the bases of the whole-variant negatives -- is curated, and every number beside
it is read from a record and the premise of the text asserted.

Workbooks are written with fixed timestamps and rewritten only when a cell
changed, so the checksum manifest stays valid across rebuilds; --check compares
cell content, not bytes, because another openpyxl version serializes the same
cells differently. The supporting texts are built first, by build_si_texts.py;
this script writes the manifest last, from the files on disk.

Usage:
    python scripts/build_supplement.py           # write the workbooks and the manifest
    python scripts/build_supplement.py --check   # rebuild in memory, compare, write nothing
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import io
import json
import pathlib
import re
import sys
from fractions import Fraction
from math import comb

import numpy as np
import openpyxl
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy.stats import fisher_exact

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"
SI_IN = REPO / "data" / "analysis_inputs" / "si"
sys.path[:0] = [str(REPO / "scripts"), str(REPO)]
import benchmark_population as bp  # noqa: E402
import build_manuscript as bm  # noqa: E402
import build_measured_effect_tables as bmet  # noqa: E402
from build_si_texts import normalize_zip  # noqa: E402

TAGS = ("t10", "t15", "t20", "t25", "tSAP")
OP_NAMES = {"t10": "1.0 kcal/mol", "t15": "1.5 kcal/mol", "t20": "2.0 kcal/mol", "t25": "2.5 kcal/mol (reference)",
            "tSAP": "Prediction-interval limit (2.9 / 2.9 / 3.5)"}
GM = {"consistent": Fraction(1), "partial": Fraction(1, 2), "inconsistent": Fraction(0)}
CLASS_NAMES = {"fold_mechanism": "Fold", "mixed_structural": "Mixed fold and interface",
               "ppi_destab_mechanism": "Interface destabilization", "ppi_stab_mechanism": "Interface stabilization",
               "structurally_silent": "Non-disruptive"}
SILENT = "structurally_silent"
EVIDENCE_LABEL = {"E1_quantitative_energetic": "E1 quantitative energetic",
                  "E2_quantitative_functional": "E2 quantitative functional",
                  "E3_qualitative_experimental": "E3 qualitative experimental",
                  "E4_population_frequency": "E4 population frequency",
                  "E5_inferred_no_axis_assay": "E5 inferred (no axis assay)"}


def J(name):
    return json.loads((RO / name).read_text())


# ------------------------------------------------------------------ rounding
def exact(x):
    """The exact rational behind a value: a float within float error of a
    fraction with a small denominator is that fraction."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, (int, np.integer)):
        return Fraction(int(x))
    f = Fraction(repr(float(x)))
    g = f.limit_denominator(10 ** 6)
    return g if abs(float(g) - float(x)) < 1e-12 else f


def fmt(x, nd):
    """Round once, half to even, on the exact value. ASCII minus, so that a
    spreadsheet reads the cell as a number."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "-"
    q = exact(x)
    neg = q < 0
    a = -q if neg else q
    scaled = a * 10 ** nd
    lo = scaled.numerator // scaled.denominator
    rem = scaled - lo
    if rem > Fraction(1, 2) or (rem == Fraction(1, 2) and lo % 2 == 1):
        lo += 1
    digits = str(lo).rjust(nd + 1, "0")
    out = digits[:-nd] + "." + digits[-nd:] if nd else digits
    return ("-" if neg and lo != 0 else "") + out


def kn(k, n, nd=3):
    return "%d/%d = %s" % (k, n, fmt(Fraction(k, n), nd))


def pfmt(p):
    """p to three decimals, or to two significant figures below 0.01."""
    if p >= 0.01:
        return fmt(p, 3)
    nd = 1
    while p * 10 ** nd < 10:
        nd += 1
    return fmt(p, nd)


def ci(pair, nd=3):
    return "%s to %s" % (fmt(pair[0], nd), fmt(pair[1], nd))


def nz(s):
    """A copied display value that rounds to zero is written without a sign."""
    return re.sub(r"^-(0\.0+)$", r"\1", s) if isinstance(s, str) else s


def vn(gene, variant):
    return "%s %s" % (str(gene).upper(), variant)


# ------------------------------------------------------------------ records
CAN = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
B = bp.benchmark(CAN).reset_index(drop=True)
assert len(B) == bp.N_BENCHMARK == 56 and len(CAN) == bp.N_SCORED
BKEYS = set(zip(B.system, B.variant))
MECH = B.expected_mech_class != SILENT
N_MECH, N_ND = int(MECH.sum()), int((~MECH).sum())
assert (N_MECH, N_ND) == (29, 27)
LED = pd.read_csv(RO / "COMAVI_evidence_ledger.csv", keep_default_na=False, dtype=str)
LB = LED[[k in BKEYS for k in zip(LED.system, LED.variant)]].reset_index(drop=True)
assert len(LED) == 99 and len(LB) == 97
CT = pd.read_csv(RO / "cohort_tables" / "COMAVI_TableS1_variants.csv", keep_default_na=False, dtype=str)
assert len(CT) == len(CAN) and not CAN.duplicated(["gene", "variant"]).any()
# keyed by (lower-case gene, variant), which is unique in the scored cohort
CTK = {(r["Gene"].lower(), r["Variant"]): r for r in CT.to_dict("records")}
CANK = {(str(r.gene).lower(), r.variant): r for _, r in CAN.iterrows()}
assert len(CTK) == len(CANK) == len(CAN)


def grades(tag, mask=None):
    g = B["mech_consistency_%s" % tag].map(GM)
    assert g.notna().all(), tag
    return g if mask is None else g[mask]


SWEEP = {r["threshold_tag"]: r for r in J("COMAVI_threshold_sweep.json")["rows"]}
EXT = {r["threshold_tag"]: r for r in J("COMAVI_threshold_sweep_extended.json")["rows"]}
for t in TAGS:
    tot = sum(grades(t))
    assert SWEEP[t]["whole_variant_n"] == 56 and Fraction(str(SWEEP[t]["whole_variant_k"])) == tot, t

# manuscript reference numbering, first appearance, as the article prints it
NUMBERS = bm.number_citations(bm.parse_source(bm.SOURCE.read_text()))
CITE_KEYS = {r["evidence_citation"]: r["reference_keys"].split(";")
             for r in csv.DictReader((SI_IN / "ledger_citation_keys.csv").open())}


def refs(keys):
    missing = [k for k in keys if k not in NUMBERS]
    assert not missing, "cited in S1 Table but not in the article: %s" % missing
    return "[%s]" % ", ".join(str(n) for n in sorted(NUMBERS[k] for k in keys))


# ================================================================ S1 Table
def s1_sheets():
    unmapped = sorted(set(LB.evidence_citation) - set(CITE_KEYS))
    assert not unmapped, "ledger citations without reference keys: %s" % unmapped
    assay = LB[LB.complex_amount_assay != ""]
    assert len(assay) == 10 and set(assay.axis) == {"binding"}
    ledger_rows = [[r.system, r.variant, r.axis, r.expected_token, r.evidence_type, r.evidence_directness,
                    r.evidence_basis, r.evidence_citation, r.complex_amount_assay or "-", refs(CITE_KEYS[r.evidence_citation])]
                   for r in LB.itertuples()]
    no_axis = sorted(k for k in BKEYS if k not in set(zip(LB.system, LB.variant)))
    wvn = list(csv.DictReader((SI_IN / "whole_variant_negatives.csv").open()))
    assert sorted((d["system"], d["variant"]) for d in wvn) == no_axis and len(no_axis) == 4
    for d in wvn:
        row = B[(B.system == d["system"]) & (B.variant == d["variant"])].iloc[0]
        assert row.expected_mech_class == SILENT
    d96 = B[B.variant == "D96V"].iloc[0]
    assert abs(float(d96.multimer_plddt_max) - 69.45) < 1e-9 and abs(float(d96.monomer_plddt) - 97.16) < 1e-9
    wvn_rows = [[d["system"], d["gene"], d["variant"], "no disruption on the modeled axes (whole-variant level)",
                 d["evidence_basis"], d["evidence_citation"], refs(d["reference_keys"].split(";"))] for d in wvn]
    return [
        dict(name="evidence_ledger",
             title="S1 Table, evidence_ledger sheet. Committed structural expectations with their evidence basis and "
                   "citation. Each row is one committed axis-level expectation of a benchmark variant; together the %d "
                   "rows are the benchmark's axis-level ground truth. complex_amount_assay names the assay for the %d "
                   "binding expectations measured only as an amount of complex. The reference column gives the number "
                   "of each cited source in the main article's reference list." % (len(LB), len(assay)),
             header=["system", "variant", "axis", "expected_token", "evidence_type", "evidence_directness",
                     "evidence_basis", "evidence_citation", "complex_amount_assay", "reference"],
             rows=ledger_rows, widths=[14, 9, 13, 9, 26, 11, 60, 48, 22, 10]),
        dict(name="whole_variant_negatives",
             title="S1 Table, whole_variant_negatives sheet. The %d non-disruptive variants graded at the "
                   "whole-variant level. Each is curated to carry no disruption on any modeled axis, but no axis-level "
                   "expectation is committed, so none appears in the evidence_ledger sheet. They are included in every "
                   "whole-variant metric; the main text and S1 Text report the score without them as a sensitivity "
                   "view." % len(wvn),
             header=["system", "gene", "variant", "expectation", "evidence_basis", "evidence_citation", "reference"],
             rows=wvn_rows, widths=[14, 8, 8, 30, 70, 20, 10]),
    ]


# ================================================================ S2 Table
MODEL_SCOPE = [
    ["Representable structural", "Known disruption maps to one or more of the energetic outputs COMAVI computes on the "
     "supplied structure.", "Graded on the complete pattern of applicable axes."],
    ["Structural context missing", "The mechanism acts through a partner, ligand, modification, oligomer or state that "
     "the supplied structure does not contain.", "Excluded from the benchmark and listed with its outputs "
     "(excluded_variants sheet). A measured axis stays in the evidence ledger but enters no reported metric."],
    ["Nonstructural or static-invisible", "Catalytic, regulatory, dynamic, aggregative or post-translational mechanism "
     "outside the static energetic outputs.", "Graded as non-disruptive when the literature supports the modeled axes, "
     "through committed neutral axes or at the whole-variant level."],
    ["Mechanistically uncommitted", "The literature commits no direction on any modeled axis.",
     "Excluded from the benchmark and listed with its outputs (excluded_variants sheet)."],
]
EXCLUSION_CLASS = {"no_committed_expectation": "Mechanistically uncommitted",
                   "mechanism_not_observable_on_structure": "Structural context missing"}

# (gene variant, axis label, measured text or None for a calibration value, interpretation, premise)
DISAGREEMENTS = [
    ("BRCA1 V1736A", "Isolated-subunit stability", None, "monomer",
     "Magnitude miss: the measured destabilization is above the reference and the prediction below it.", "miss"),
    ("MLH1 R755W", "Isolated-subunit stability", "Curated destabilizing", "monomer",
     "Miss against the curated fold mechanism; no stability axis reaches the threshold.", "curated"),
    ("BRCA1 R1751Q", "Isolated-subunit stability", None, "monomer",
     "Both values are below the reference; the expectation is destabilizing because the measured change is positive.",
     "both_below"),
    ("BRCA1 L1664P", "Isolated-subunit stability", None, "monomer",
     "Both values are below the reference; the expectation is destabilizing because the measured change is positive.",
     "both_below"),
    ("HBB N102T", "Binding", None, "binding",
     "Magnitude miss: the measured destabilization is above the reference and the predicted binding change is near "
     "zero.", "miss"),
    ("HBB W37Y", "Binding", None, "binding",
     "Both values are below the reference; the expectation is destabilizing because the measured change is positive.",
     "both_below"),
    ("MSH2 A636P", "Isolated-subunit and assembled-complex stability", "Neutral on both modeled fold axes", "both_fold",
     "False fold calls against literature-supported neutral expectations; the mismatch-processing defect lies outside "
     "the modeled axes.", "false_fold"),
]


def _pred(key, col):
    v = CTK[key][col]
    return None if v in ("", "-") else float(v)


def exact_perm_p(system, truth, strong):
    """Exact within-system permutation p (convolved hypergeometrics), as a Fraction."""
    pop = pd.DataFrame({"system": list(system), "truth": list(truth), "strong": list(strong)})
    obs = int((pop["strong"] & pop["truth"]).sum())
    tot = {0: Fraction(1)}
    for _, g in pop.groupby("system"):
        n_, k_, s_ = len(g), int(g["truth"].sum()), int(g["strong"].sum())
        d = {x: Fraction(comb(s_, x) * comb(n_ - s_, k_ - x), comb(n_, k_))
             for x in range(max(0, k_ - (n_ - s_)), min(k_, s_) + 1)}
        nt = {}
        for a, pa in tot.items():
            for x, px in d.items():
                nt[a + x] = nt.get(a + x, Fraction(0)) + pa * px
        tot = nt
    return sum(p for s, p in tot.items() if s >= obs)


def s2_sheets():
    # ---- benchmark_variants
    bv = []
    stable_change = {"four": 0, "five": 0}
    for _, r in B.iterrows():
        key = (str(r.gene).lower(), r.variant)
        c = CTK[key]
        assert c["COMAVI agreement (t=2.5)"] == r.mech_consistency_t25
        codes = c["Agreement across thresholds"].split("/")
        want = ["C" if r["mech_consistency_%s" % t] == "consistent" else "P" if r["mech_consistency_%s" % t] == "partial"
                else "X" for t in TAGS]
        assert codes == want, (key, codes, want)
        stable_change["four"] += len(set(codes[:4])) > 1
        stable_change["five"] += len(set(codes)) > 1
        am = r.get("AM pathogenicity")
        bv.append([c["System"], str(r.gene).upper(), r.variant, c["Clinical annotation"], c["Expected mechanism"],
                   c["Structural tier"], r.mech_t25, r.mech_consistency_t25, "/".join(codes),
                   "yes" if len(set(codes)) == 1 else "no", c["Committed axes"], c["Evidence class"], c["Evidence basis"],
                   c["Citation"], nz(c["Predicted DDG, isolated subunit (kcal/mol)"]), nz(c["Predicted DDG, in complex (kcal/mol)"]),
                   nz(c["Predicted DDG, binding (kcal/mol)"]), nz(c["Measured DDG (kcal/mol)"]), c["Measured axis"], c["Structure"],
                   "n/a" if pd.isna(am) else fmt(am, 3)])
    assert (stable_change["four"], stable_change["five"]) == (14, 17), stable_change
    # ---- excluded_variants
    ex = bp.excluded(CAN)
    assert len(ex) == 5
    exr = []
    for _, r in ex.iterrows():
        key = (str(r.gene).lower(), r.variant)
        c = CTK[key]
        assert pd.isna(r.mech_consistency_t25)
        exr.append([c["System"], str(r.gene).upper(), r.variant, EXCLUSION_CLASS[r.exclusion_reason_code],
                    r.exclusion_reason, c["Clinical annotation"], r.mech_t25,
                    nz(c["Predicted DDG, isolated subunit (kcal/mol)"]), nz(c["Predicted DDG, in complex (kcal/mol)"]),
                    nz(c["Predicted DDG, binding (kcal/mol)"]), c["Structure"]])
    # ---- evidence_composition
    comp = []
    notes = {
        "E1_quantitative_energetic": "Quantitative energetic measurement on the axis, such as a melting temperature, "
                                     "an unfolding free energy or a dissociation constant.",
        "E2_quantitative_functional": "Quantitative functional readout coupled to the axis.",
        "E3_qualitative_experimental": "Qualitative experimental evidence on the axis.",
        "E4_population_frequency": "Population frequency incompatible with the disruption.",
        "E5_inferred_no_axis_assay": "Axis state inferred from a mechanistic or retained-function result without an "
                                     "axis-specific assay; never assigned from missing data alone.",
    }
    for et, note in notes.items():
        sub = LB[LB.evidence_type == et]
        tok = sub.expected_token.value_counts().to_dict()
        assert set(tok) <= {"neutral", "destab", "stab"}
        d = sub.evidence_directness.value_counts().to_dict()
        extra = ""
        if tok.get("stab"):
            gene_of = dict(zip(zip(CAN.system, CAN.variant), CAN.gene.astype(str).str.upper()))
            axis_name = {"binding": "binding", "monomer": "isolated-subunit stability",
                         "fold_complex": "assembled-complex stability"}
            st = sub[sub.expected_token == "stab"]
            extra = " The stabilizing expectation is %s." % "; ".join(
                "%s %s %s" % (gene_of[(s, v)], v, axis_name[a]) for s, v, a in zip(st.system, st.variant, st.axis))
        if not tok.get("destab") and not tok.get("stab"):
            extra = " Every expectation in this class is neutral."
        comp.append([EVIDENCE_LABEL[et], len(sub), tok.get("neutral", 0),
                     tok.get("destab", 0), tok.get("stab", 0), d.get("direct", 0), d.get("coupled", 0),
                     d.get("inferred", 0), note + extra])
    assert sum(r[1] for r in comp) == len(LB)
    comp.append(["All classes", len(LB), int((LB.expected_token == "neutral").sum()), int((LB.expected_token == "destab").sum()),
                 int((LB.expected_token == "stab").sum()), int((LB.evidence_directness == "direct").sum()),
                 int((LB.evidence_directness == "coupled").sum()), int((LB.evidence_directness == "inferred").sum()), ""])
    # ---- context_comparators
    cbt = pd.read_csv(RO / "isds_v1" / "ISDS_v1_component_binary_baselines.csv").set_index("screen")
    pvv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
    assert len(pvv) == 46
    truth = pvv.structural_ground_truth.astype(bool)
    screens = [("interface_status_alone", "Interface status alone (any confidence-passing interface contact)",
                pvv.interface_only_positive.astype(bool)),
               ("tier_without_interface_bonus", "Tier 1-2 without the interface bonus", pvv.no_interface_score >= 3.0),
               ("full_tier_1_2", "Tier 1-2 (full tier)", pvv.reported_tier.isin(["Tier 1", "Tier 2"]))]
    cc = []
    for key, name, strong in screens:
        q = cbt.loc[key]
        tp, fn, fp, tn = int(q.tp), int(q.fn), int(q.fp), int(q.tn)
        assert (tp, fp) == (int((strong & truth).sum()), int((strong & ~truth).sum())), key
        assert tp + fn == int(truth.sum()) and fp + tn == int((~truth).sum())
        pp = exact_perm_p(pvv.system, truth, strong)
        assert abs(float(pp) - float(q.within_system_permutation_p)) < 5e-5, key
        fp_ = fisher_exact([[tp, fn], [fp, tn]])[1]
        # p to three decimals, as in S2 Text Table E, which shows the same screens
        cc.append([name, "%d/%d" % (tp, fn), "%d/%d" % (fp, tn), fmt(Fraction(tp, tp + fn), 3), fmt(Fraction(tn, fp + tn), 3),
                   fmt(Fraction(tp, tp + fp), 3), fmt((Fraction(tp, tp + fn) + Fraction(tn, fp + tn)) / 2, 3),
                   fmt(fp_, 3), fmt(pp, 3)])
    # ---- fold_binding_disagreements
    pts = pd.read_csv(RO / "COMAVI_delta_calibration_points.csv")
    measured = dict(zip(pts.variant, pts.measured_kcal))
    fb = []
    for name, axis_label, mtext, axis, interp, premise in DISAGREEMENTS:
        gene, var = name.split()
        key = (gene.lower(), var)
        r = CANK[key]
        assert (r.system, var) in BKEYS
        if axis == "both_fold":
            pm, pc = _pred(key, "Predicted DDG, isolated subunit (kcal/mol)"), _pred(key, "Predicted DDG, in complex (kcal/mol)")
            led = LB[(LB.system == r.system) & (LB.variant == var)].set_index("axis").expected_token
            assert led["monomer"] == "neutral" and led["fold_complex"] == "neutral" and pm >= 2.5 and pc >= 2.5
            fb.append([name, axis_label, mtext, "%s / %s" % (fmt(pm, 2), fmt(pc, 2)), interp])
            continue
        col = {"monomer": "Predicted DDG, isolated subunit (kcal/mol)", "binding": "Predicted DDG, binding (kcal/mol)"}[axis]
        pred = _pred(key, col)
        led = LB[(LB.system == r.system) & (LB.variant == var)].set_index("axis").expected_token
        assert led[axis] == "destab", name
        if premise == "curated":
            assert pred < 2.5
            fb.append([name, axis_label, mtext, fmt(pred, 2), interp])
            continue
        meas = float(measured[var])
        if premise == "miss":
            assert meas >= 2.5 > pred, name
        elif premise == "both_below":
            assert 0 < meas < 2.5 and pred < 2.5, name
        fb.append([name, axis_label, fmt(meas, 2), fmt(pred, 2), interp])
    return [
        dict(name="benchmark_variants",
             title="S2 Table, benchmark_variants sheet. The %d benchmark variants. Each row gives the clinical "
                   "annotation, the expected mechanism defined from primary literature, the structural-context tier, the "
                   "COMAVI call and grade at the 2.5 kcal/mol reference, the grade at every operating point (C "
                   "consistent, P partial, X inconsistent, at 1.0, 1.5, 2.0 and 2.5 kcal/mol and at the per-axis "
                   "prediction-interval limit), the committed axes with their evidence class, basis and citation, the "
                   "predicted and measured \u0394\u0394G (kcal/mol), the structure used and the AlphaMissense "
                   "pathogenicity score where AlphaMissense covers the substitution. %d variants change grade across the "
                   "four uniform thresholds and %d across all five operating points." % (len(bv), stable_change["four"],
                                                                                         stable_change["five"]),
             header=["System", "Gene", "Variant", "Clinical annotation", "Expected mechanism", "Structural tier",
                     "COMAVI call at 2.5 kcal/mol", "Grade at 2.5 kcal/mol", "Grades at the five operating points",
                     "Same grade at every operating point", "Committed axes", "Evidence class", "Evidence basis",
                     "Citation", "Predicted \u0394\u0394G, isolated subunit", "Predicted \u0394\u0394G, assembled complex",
                     "Predicted \u0394\u0394G, binding", "Measured \u0394\u0394G", "Measured axis", "Structure",
                     "AlphaMissense pathogenicity"],
             rows=bv, widths=[30, 8, 9, 14, 24, 9, 34, 12, 12, 10, 24, 10, 60, 40, 10, 10, 10, 9, 14, 16, 12]),
        dict(name="evidence_composition",
             title="S2 Table, evidence_composition sheet. The %d committed expectations of the benchmark variants by "
                   "evidence class, expected direction and directness of the evidence." % len(LB),
             header=["Evidence class", "Expectations", "Neutral", "Destabilizing", "Stabilizing", "Direct", "Coupled",
                     "Inferred", "Definition"],
             rows=comp, widths=[28, 12, 9, 12, 10, 8, 8, 9, 80]),
        dict(name="model_scope",
             title="S2 Table, model_scope sheet. Model-scope rules: how a variant is treated when the supplied structure "
                   "does or does not contain the context its known mechanism needs.",
             header=["Class", "Definition", "Treatment"], rows=MODEL_SCOPE, widths=[30, 70, 70]),
        dict(name="excluded_variants",
             title="S2 Table, excluded_variants sheet. The %d scored candidates excluded from the benchmark because no "
                   "gradeable expectation could be defined for them, with the reason and the COMAVI outputs at the 2.5 "
                   "kcal/mol reference (\u0394\u0394G in kcal/mol). They enter no reported metric." % len(exr),
             header=["System", "Gene", "Variant", "Model-scope class", "Reason", "Clinical annotation",
                     "COMAVI call at 2.5 kcal/mol", "Predicted \u0394\u0394G, isolated subunit",
                     "Predicted \u0394\u0394G, assembled complex", "Predicted \u0394\u0394G, binding", "Structure"],
             rows=exr, widths=[30, 8, 9, 26, 70, 14, 30, 10, 10, 10, 16]),
        dict(name="context_comparators",
             title="S2 Table, context_comparators sheet. Structural-context screens on the %d-variant prioritization "
                   "population (%d structural mechanisms, %d non-disruptive variants, %d systems). Within-system "
                   "permutation p values are exact." % (len(pvv), int(truth.sum()), int((~truth).sum()), pvv.system.nunique()),
             header=["Screen", "Mechanisms recovered / missed", "Non-disruptive called / rejected", "Sensitivity",
                     "Specificity", "Positive predictive value", "Balanced accuracy", "Fisher exact p",
                     "Within-system permutation p"],
             rows=cc, widths=[44, 14, 14, 11, 11, 12, 11, 11, 13]),
        dict(name="fold_binding_disagreements",
             title="S2 Table, fold_binding_disagreements sheet. Directly measured or curated fold and binding "
                   "disagreements at the 2.5 kcal/mol reference (\u0394\u0394G in kcal/mol; measured values are the "
                   "calibration points of S2 Fig).",
             header=["Variant", "Axis", "Measured or curated expectation", "Predicted \u0394\u0394G", "Interpretation"],
             rows=fb, widths=[14, 30, 26, 14, 90]),
    ]


# ================================================================ S3 Table
def s3_sheets():
    # ---- threshold_sweep
    cols = [("Detection", "detection_k", "detection_n"), ("Attribution", "attribution_k", "attribution_n"),
            ("Correct rejection", "correct_rejection_k", "correct_rejection_n"),
            ("Whole-variant score", "whole_variant_k", "whole_variant_n"),
            ("Axis: isolated subunit", "axis_monomer_k", "axis_monomer_n"),
            ("Axis: assembled complex", "axis_fold_k", "axis_fold_n"),
            ("Axis: binding", "axis_binding_k", "axis_binding_n"),
            ("Axis: either assembled axis (amount of complex)", "axis_assembled_either_k", "axis_assembled_either_n"),
            ("Structural agreement", "structural_agreement_k", "structural_agreement_n"),
            ("Detection, isolated subunit only", "detection_isolated_k", "detection_n"),
            ("Correct rejection, isolated subunit only", "correct_rejection_isolated_k", "correct_rejection_n")]
    sweep_rows = []
    for t in TAGS + ("t30", "t35", "t40"):
        rec = SWEEP.get(t) or EXT[t]
        if t in SWEEP:
            for k in rec:
                if k in EXT[t] and not isinstance(rec[k], (dict, list)) and rec[k] == rec[k]:
                    assert EXT[t][k] == rec[k], (t, k)
        label = {"t25": "2.5 (reference)", "tSAP": "2.9/2.9/3.5"}.get(t, str(rec["threshold_kcal_mol"]))
        if t in ("t30", "t35", "t40"):
            label += " (post hoc)"
        row = [label]
        for _, kk, nn in cols:
            k, n = rec.get(kk), rec.get(nn)
            if k is None or n is None or (isinstance(k, float) and np.isnan(k)):
                row.append("-")
            else:
                row.append(fmt(Fraction(str(k)) / int(n), 3))
        sweep_rows.append(row)
    t25 = SWEEP["t25"]
    denoms = "; ".join("%s n = %d" % (h.lower(), t25[nn]) for h, _, nn in cols[:9])
    # ---- operating_point_sweep
    opd = J("COMAVI_operating_point_diagnostics.json")
    me = J("COMAVI_measured_effect_tables.json")
    rec_by = {d["threshold"]: d for d in me["measured_effects_recovered"]["by_threshold"]}
    n_meas = me["measured_effects_recovered"]["population"]["n"]
    assert bmet.DESTABILIZER_MIN == 1.0 and me["measured_effects_recovered"]["population"]["definition"] == \
        "measured >= +1.0 kcal/mol; recovered when predicted >= +threshold (same sign)"
    op = []
    for t in TAGS:
        am, asl = Fraction(sum(grades(t, MECH))) / N_MECH, Fraction(sum(grades(t, ~MECH))) / N_ND
        pe = opd["point_estimates"][t]
        assert abs(float(am) - pe["recovery"]) < 6e-5 and abs(float(asl) - pe["correct_rejection"]) < 6e-5
        pb = opd["paired_bootstrap"]
        pooled = Fraction(sum(grades(t))) / 56
        mr = rec_by[{"t10": "1.0", "t15": "1.5", "t20": "2.0", "t25": "2.5", "tSAP": "tSAP"}[t]]
        assert mr["n"] == n_meas
        op.append([OP_NAMES[t], "%s (%s)" % (fmt(am, 3), ci(pb["recovery"][t]["ci95"])),
                   "%s (%s)" % (fmt(asl, 3), ci(pb["correct_rejection"][t]["ci95"])),
                   "%s/56 = %s" % (fmt(sum(grades(t)), 1), fmt(pooled, 3)), fmt((am + asl) / 2, 3),
                   "%d/%d" % (mr["recovered"], mr["n"])])
    # ---- missed_mechanisms
    dva = J("COMAVI_detection_vs_attribution.json")
    msd = dva["missed_at_reference"]
    assert msd["n_missed"] == len(msd["variants"]) == 9 and msd["threshold"] in (2.5, "2.5")
    mm = []
    for d in msd["variants"]:
        meas = d.get("measured_kcal_mol")
        rel = "-" if meas is None else ("below the reference" if meas < 2.5 else "above the reference")
        r = CANK[(d["gene"].lower(), d["variant"])]
        # no axis fires: the call is either silence or an energy-neutral context label
        assert (r.mech_t25 == "No structural effect detected" or "(DDG neutral)" in r.mech_t25) \
            and r.expected_mech_class == d["expected_class"], (d["variant"], r.mech_t25)
        mm.append([d["gene"], d["variant"], CLASS_NAMES[d["expected_class"]], r.mech_t25, fmt(d["predicted_monomer"], 2),
                   fmt(d["predicted_complex_fold"], 2), fmt(d["predicted_binding"], 2),
                   "-" if meas is None else fmt(meas, 2), rel,
                   "-" if meas is None else re.sub(r"^evidence ledger \((\w+); (.*)\)$", r"\2 (\1)", d["measured_source"])])
    below = sorted("%s %s" % (d["gene"], d["variant"]) for d in msd["variants"]
                   if d.get("measured_kcal_mol") is not None and d["measured_kcal_mol"] < 2.5)
    assert below == sorted(msd["measured_below_reference"]) and len(below) == 4
    # ---- disagreement_classes
    cr = me["call_relationships"]
    win = bmet.BORDERLINE_WINDOW
    defs = {"Agreement": "Measured and predicted calls agree at 2.5 kcal/mol.",
            "Borderline disagreement": "Calls differ at 2.5 kcal/mol, and both values lie within %s kcal/mol of it." % fmt(win, 1),
            "Substantive disagreement": "Calls differ at 2.5 kcal/mol, and at least one value lies more than %s kcal/mol "
                                        "from it." % fmt(win, 1)}
    wo = {c["call_relationship"]: c for c in cr["sensitivity_excluding_glu73"]["classes"]}
    dc = [[c["call_relationship"], c["n"], fmt(c["median_abs_error"], 2), wo[c["call_relationship"]]["n"],
           fmt(wo[c["call_relationship"]]["median_abs_error"], 2), defs[c["call_relationship"]], "; ".join(c["examples"])]
          for c in cr["classes"]]
    assert sum(r[1] for r in dc) == cr["n_comparisons"] == 63
    assert sum(r[3] for r in dc) == cr["sensitivity_excluding_glu73"]["n_comparisons"] == 56
    # ---- tier_distribution
    pvv = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
    tn = pvv.reported_tier.astype(str).str.extract(r"(\d)").astype(float)[0]
    y = pvv.structural_ground_truth.astype(bool)
    td = []
    for label, mask in (("Tier 1", tn == 1), ("Tier 2", tn == 2), ("Tier 3", tn == 3), ("Tier 4", tn == 4),
                        ("Tiers 3-4", tn.isin([3, 4]))):
        k, n = int((mask & y).sum()), int(mask.sum())
        td.append([label, k, n - k, kn(k, n)])
    assert sum(r[1] + r[2] for r in td[:4]) == len(pvv) == 46
    # ---- evidence_states and detection_metrics
    est = J("isds_v1/ISDS_v1_evidence_state_tables.json")
    assert est["population"] == {"n": 46, "structural": 21, "no_lesion": 25}
    es, dm = [], []
    for row in est["four_evidence_states"]:
        cells = [row[s] for s in ("Convergent", "Energy only", "Context only", "Neither")]
        assert sum(c["structural"] for c in cells) == 21 and sum(c["no_lesion"] for c in cells) == 25
        es.append([OP_NAMES[row["threshold"]]] + ["%d/%d" % (c["structural"], c["no_lesion"]) for c in cells])
    for row in est["detection_metrics"]:
        out = [OP_NAMES[row["threshold"]]]
        for key in ("any_energetic_firing", "convergent_energy_and_tier"):
            q = row[key]
            assert q["tp"] + q["fn"] == 21 and q["fp"] + q["tn"] == 25
            out += ["%d/%d/%d/%d" % (q["tp"], q["fn"], q["fp"], q["tn"]),
                    "%s / %s / %s" % (fmt(Fraction(q["tp"], q["tp"] + q["fn"]), 3), fmt(Fraction(q["tn"], q["tn"] + q["fp"]), 3),
                                      fmt(Fraction(q["tp"], q["tp"] + q["fp"]), 3))]
        dm.append(out)
    # ---- pathogenicity_ranking
    pr = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_pathogenicity_comparison.csv")
    common = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_common_set.csv")
    from sklearn.metrics import roc_auc_score
    interp = {"AM pathogenicity": ("AlphaMissense", "Pathogenicity ranking"),
              "isds_v1": ("Priority score", "Structural-disruption priority"),
              "isds_energy_component": ("Energy component of the priority score", "Energetic structural evidence"),
              "isds_context_component": ("Context component of the priority score", "Structural-context strength")}
    prr = []
    for _, r in pr.iterrows():
        assert r.n == len(common) and r.n_pathogenic == int(common.clinical_y.sum())
        auc = roc_auc_score(common.clinical_y, common[r.descriptor])
        npairs = int(common.clinical_y.sum() * (1 - common.clinical_y).sum())
        assert abs(auc * npairs * 2 - round(auc * npairs * 2)) < 1e-9 and abs(auc - r.pathogenicity_auc) < 1e-12
        name, use = interp[r.descriptor]
        prr.append([name, "%d (%d / %d)" % (r.n, r.n_pathogenic, r.n - r.n_pathogenic),
                    fmt(Fraction(int(round(auc * npairs * 2)), 2 * npairs), 3), use])
    # ---- robustness
    z = np.load(RO / "stress_tests" / "comavi_stress_draws.npz")
    vb = J("stress_tests/COMAVI_stress_verified_statistics.json")
    stc = pd.read_csv(RO / "stress_tests" / "comavi_stress_tests.csv").set_index("test")
    sa_k, sa_n = int(B.structural_agreement_n_t25.sum()), int(B.structural_agreement_d_t25.sum())
    assert (sa_k, sa_n) == tuple(vb["population"]["primary_structural_agreement"]) == (t25["structural_agreement_k"],
                                                                                        t25["structural_agreement_n"])
    mc_obs, sa_obs = Fraction(sum(grades("t25"))) / 56, Fraction(sa_k, sa_n)
    assert abs(float(mc_obs) - float(z["mc_obs"][0])) < 1e-12 and abs(float(sa_obs) - float(z["sa_obs"][0])) < 1e-9
    g25 = B.assign(g=grades("t25"), sk=B.structural_agreement_n_t25, sn=B.structural_agreement_d_t25)
    per = g25.groupby("system").agg(g=("g", "sum"), n=("g", "size"), sk=("sk", "sum"), sn=("sn", "sum"))
    loso_mc = [(sum(grades("t25")) - Fraction(per.loc[s, "g"])) / (56 - int(per.loc[s, "n"])) for s in per.index]
    loso_sa = [Fraction(sa_k - int(per.loc[s, "sk"]), sa_n - int(per.loc[s, "sn"])) for s in per.index]
    lcsv = pd.read_csv(RO / "stress_tests" / "comavi_leave_one_system_out.csv")
    assert abs(float(min(loso_mc)) - lcsv.mc.min()) < 6e-5 and abs(float(max(loso_mc)) - lcsv.mc.max()) < 6e-5
    assert abs(float(min(loso_sa)) - lcsv.sa_val.min()) < 6e-5 and abs(float(max(loso_sa)) - lcsv.sa_val.max()) < 6e-5

    def null_txt(arr, p):
        # p to four decimals, as the main text and S2 Text report it
        return "%s (SD %s); p = %s" % (fmt(float(arr.mean()), 3), fmt(float(arr.std()), 3), fmt(p, 4))
    p_mc = Fraction(int((z["mc_null"] >= float(mc_obs) - 1e-12).sum()) + 1, len(z["mc_null"]) + 1)
    p_sa = Fraction(int((z["sa_null"] >= float(sa_obs) - 1e-12).sum()) + 1, len(z["sa_null"]) + 1)
    assert abs(float(p_mc) - stc.loc["permutation null (MC)", "p_value"]) < 5e-5
    assert abs(float(p_sa) - stc.loc["permutation null (SA)", "p_value"]) < 5e-5
    rb = [["Observed at 2.5 kcal/mol", "%s/56 = %s" % (fmt(sum(grades("t25")), 1), fmt(mc_obs, 3)), kn(sa_k, sa_n)],
          ["Label-permutation null, %s permutations: mean (SD); p" % format(len(z["mc_null"]), ","),
           null_txt(z["mc_null"], float(p_mc)), null_txt(z["sa_null"], float(p_sa))],
          ["Leave one system out (%d systems): range" % len(per), "%s to %s" % (fmt(min(loso_mc), 3), fmt(max(loso_mc), 3)),
           "%s to %s" % (fmt(min(loso_sa), 3), fmt(max(loso_sa), 3))],
          ["Replicate noise, %d draws: mean (SD)" % len(z["noise_mc"]),
           "%s (SD %s)" % (fmt(float(z["noise_mc"].mean()), 3), fmt(float(z["noise_mc"].std()), 3)),
           "%s (SD %s)" % (fmt(float(z["noise_sa"].mean()), 3), fmt(float(z["noise_sa"].std()), 3))],
          ["Variant bootstrap, %s draws: 95%% interval" % format(vb["draws"]["variant_bootstrap"], ","),
           ci(vb["variant_bootstrap_ci95"]["mc"]["full_precision"]), ci(vb["variant_bootstrap_ci95"]["sa"]["full_precision"])],
          ["System-cluster bootstrap, %s draws: 95%% interval" % format(vb["draws"]["system_cluster_bootstrap"], ","),
           ci(vb["system_cluster_ci95"]["mc"]["full_precision"]), ci(vb["system_cluster_ci95"]["sa"]["full_precision"])]]
    # ---- skempi_validation
    skv = J("COMAVI_skempi_validation.json")
    cx = skv["complexes"]
    bbv = pd.read_csv(REPO / "supplement" / "skempi" / "bb_binding_validation.csv")
    jtv = pd.read_csv(REPO / "supplement" / "skempi" / "jt_binding_validation.csv")
    g73 = bbv.foldx_code.str.match(r"E[A-Z]?73")  # FoldX numbering, as build_skempi_validation.GLU73
    assert int(g73.sum()) == skv["glu73_cluster"]["n"] == 7

    def meth(df):
        return "; ".join("%s %d" % (m, c) for m, c in sorted(df.method.value_counts().items(), key=lambda kv: (-kv[1], kv[0])))
    assert len(bbv) == cx["barnase_barstar_all"]["n"] and len(jtv) == cx["tem1_blip"]["n"]
    sk = []
    for key, name, df in (("barnase_barstar_all", "Barnase-barstar, all", bbv),
                          ("barnase_barstar_excluding_glu73", "Barnase-barstar, without Glu73", bbv[~g73]),
                          ("tem1_blip", "TEM1-BLIP, BLIP side", jtv)):
        d = cx[key]
        assert len(df) == d["n"]
        sk.append([name, meth(df), d["n"], "%s (p = %s)" % (fmt(d["pearson_r"], 2), pfmt(d["pearson_p"])),
                   "%s (p = %s)" % (fmt(d["spearman_rho"], 2), pfmt(d["spearman_p"])),
                   "%d/%d" % (d["sign_correct"], d["destabilizing_n"])])
    for key, name in (("combined_including_glu73", "Both complexes"), ("combined_excluding_glu73", "Both complexes, without Glu73")):
        d = skv[key]
        sk.append([name, "mixed", d["n"], "%s (p = %s)" % (fmt(d["pearson_r"], 2), pfmt(d["pearson_p"])),
                   "%s (p = %s)" % (fmt(d["spearman_rho"], 2), pfmt(d["spearman_p"])),
                   "%d/%d" % (d["sign_correct"], d["destabilizing_n"])])
    return [
        dict(name="threshold_sweep",
             title="S3 Table, threshold_sweep sheet. Every graded metric at every decision threshold. The first five "
                   "rows are the prespecified operating points; the three post hoc rows (3.0, 3.5 and 4.0 kcal/mol) "
                   "were evaluated after those results were known, and the whole-variant score is not computed at them. "
                   "Denominators at 2.5 kcal/mol: %s. The either-assembled-axis column holds the expectations of a "
                   "change measured only as an amount of complex, scored jointly." % denoms,
             header=["Threshold (kcal/mol)"] + [h for h, _, _ in cols], rows=sweep_rows,
             widths=[22] + [12] * len(cols)),
        dict(name="operating_point_sweep",
             title="S3 Table, operating_point_sweep sheet. Complete operating-point sweep. The mechanism arm is the mean "
                   "grade over the %d structural-mechanism variants and the non-disruptive arm over the %d non-disruptive "
                   "variants, with 95%% intervals from %s paired resamples of whole systems; the pooled score is the "
                   "mechanism-pattern score over all 56; the balanced score is the mean of the two arms. A measured "
                   "destabilization (measured \u0394\u0394G of at least 1.0 kcal/mol, n = %d) is recovered when the "
                   "prediction reaches the threshold with the same sign." % (N_MECH, N_ND,
                                                                            format(opd["paired_bootstrap"]["n_draws"], ","), n_meas),
             header=["Operating point", "Mechanism arm (n = %d)" % N_MECH, "Non-disruptive arm (n = %d)" % N_ND,
                     "Pooled score (n = 56)", "Balanced score", "Measured destabilizations recovered (n = %d)" % n_meas],
             rows=op, widths=[40, 22, 22, 18, 12, 18]),
        dict(name="missed_mechanisms",
             title="S3 Table, missed_mechanisms sheet. The %d structural mechanisms on which no axis reaches the 2.5 "
                   "kcal/mol reference, with the predicted \u0394\u0394G on each axis and, where one exists, the measured "
                   "\u0394\u0394G (kcal/mol). For %d of them the measured effect is itself below the reference." % (len(mm), len(below)),
             header=["Gene", "Variant", "Expected mechanism", "COMAVI call at 2.5 kcal/mol", "Predicted \u0394\u0394G, isolated subunit",
                     "Predicted \u0394\u0394G, assembled complex", "Predicted \u0394\u0394G, binding", "Measured \u0394\u0394G",
                     "Measured effect", "Measurement source"],
             rows=mm, widths=[8, 9, 26, 34, 12, 12, 12, 10, 20, 60]),
        dict(name="disagreement_classes",
             title="S3 Table, disagreement_classes sheet. Call relationships between measured and predicted energies at "
                   "the 2.5 kcal/mol reference, for all %d comparisons and without the seven barnase Glu73 substitutions."
                   % cr["n_comparisons"],
             header=["Call relationship", "n", "Median absolute error (kcal/mol)", "n without Glu73",
                     "Median without Glu73 (kcal/mol)", "Definition", "Examples"],
             rows=dc, widths=[24, 6, 14, 10, 14, 70, 26]),
        dict(name="tier_distribution",
             title="S3 Table, tier_distribution sheet. Structural mechanisms by structural-context tier on the 46-variant "
                   "prioritization population.",
             header=["Tier", "Structural mechanism", "Non-disruptive", "Structural fraction"], rows=td, widths=[12, 12, 14, 20]),
        dict(name="evidence_states",
             title="S3 Table, evidence_states sheet. Four evidence states at each operating point on the 46-variant "
                   "prioritization population. Each cell is structural mechanisms / non-disruptive variants; convergent "
                   "is an energy call with Tier 1-2 context.",
             header=["Operating point", "Convergent", "Energy only", "Context only", "Neither"], rows=es,
             widths=[40, 12, 12, 12, 12]),
        dict(name="detection_metrics",
             title="S3 Table, detection_metrics sheet. Detection of structural mechanisms on the 46-variant "
                   "prioritization population by any energy call and by a convergent call (energy and Tier 1-2). Counts "
                   "are TP/FN/FP/TN; rates are sensitivity / specificity / positive predictive value.",
             header=["Operating point", "Any energy call: counts", "Any energy call: rates", "Convergent: counts",
                     "Convergent: rates"], rows=dm, widths=[40, 14, 24, 14, 24]),
        dict(name="pathogenicity_ranking",
             title="S3 Table, pathogenicity_ranking sheet. Descriptive ranking of pathogenic against benign variants on "
                   "the %d prioritization variants that carry a clinical label and an AlphaMissense score. The priority "
                   "score is not a pathogenicity predictor; the comparison is descriptive." % len(common),
             header=["Score", "Variants (pathogenic / benign)", "ROC AUC for pathogenicity", "What the score measures"],
             rows=prr, widths=[40, 18, 14, 36]),
        dict(name="robustness",
             title="S3 Table, robustness sheet. Robustness of the mechanism-pattern score and the four-output structural "
                   "agreement at the 2.5 kcal/mol reference.",
             header=["Test", "Mechanism-pattern score", "Structural agreement"], rows=rb, widths=[52, 30, 30]),
        dict(name="skempi_validation",
             title="S3 Table, skempi_validation sheet. FoldX binding \u0394\u0394G against SKEMPI 2.0 measurements for "
                   "every substitution and without the barnase Glu73 cluster. Methods are the SKEMPI 2.0 method codes "
                   "with counts (ITC, isothermal titration calorimetry; SFFL, stopped-flow fluorescence; IAFL, "
                   "fluorescence inhibition assay; IASP, spectroscopic inhibition assay; SPR, surface plasmon "
                   "resonance). Barnase residues are numbered as in the structure scored by FoldX, in which the "
                   "catalytic glutamate is Glu73.",
             header=["Complex", "Methods", "Substitutions", "Pearson r", "Spearman \u03c1",
                     "Sign correct on measured destabilizations"], rows=sk, widths=[30, 30, 12, 18, 18, 16]),
    ]


# ================================================================ writing
TITLES = {"S1": "Committed structural expectations with their evidence basis and citation.",
          "S2": "Benchmark cohort, evidence and comparator tables.",
          "S3": "Operating points, evidence states and robustness tables."}


def si_captions():
    """The manuscript's @si block: item -> (title, legend)."""
    out = {}
    for line in bm.SOURCE.read_text().splitlines():
        if line.startswith("@si "):
            item, title, legend = [x.strip() for x in line[4:].split("|", 2)]
            out[item] = (title, legend)
    return out


def workbook_bytes(item, sheets):
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    idx = wb.create_sheet("index")
    head = ["Supporting information item", "sheet", "caption"]
    for j, h in enumerate(head, 1):
        c = idx.cell(1, j, h)
        c.font = Font(bold=True)
        c.fill = PatternFill("solid", fgColor="D9D9D9")
    for i, sh in enumerate(sheets, 2):
        idx.cell(i, 1, item)
        idx.cell(i, 2, sh["name"])
        idx.cell(i, 3, sh["title"]).alignment = Alignment(wrap_text=True, vertical="top")
    for col, w in zip("ABC", (24, 26, 120)):
        idx.column_dimensions[col].width = w
    for sh in sheets:
        ws = wb.create_sheet(sh["name"])
        ws.cell(1, 1, sh["title"]).font = Font(bold=True)
        for j, h in enumerate(sh["header"], 1):
            c = ws.cell(3, j, h)
            c.font = Font(bold=True)
            c.fill = PatternFill("solid", fgColor="D9D9D9")
            c.alignment = Alignment(wrap_text=True, vertical="top")
        for i, row in enumerate(sh["rows"], 4):
            assert len(row) == len(sh["header"]), (sh["name"], row)
            for j, v in enumerate(row, 1):
                ws.cell(i, j, v.item() if isinstance(v, np.generic) else v)
        for j, w in enumerate(sh.get("widths") or [], 1):
            ws.column_dimensions[get_column_letter(j)].width = w
        ws.freeze_panes = "A4"
    wb.properties.creator = "COMAVI build_supplement.py"
    buf = io.BytesIO()
    wb.save(buf)
    return normalize_zip(buf.getvalue())


def si_dir():
    hits = sorted(glob.glob(str(REPO / "submission" / "COMAVI_*_SI")))
    if len(hits) != 1:
        raise SystemExit("FAIL: expected exactly one submission/COMAVI_*_SI directory, found %s" % hits)
    return pathlib.Path(hits[0])


def body_text():
    t = (REPO / "docs" / "COMAVI_manuscript_current.txt").read_text()
    cut = [i for i in (t.find("\nSupporting information"), t.find("[SI S1_Text]")) if i > 0]
    return t[:min(cut)] if cut else t


def manifest(SI, books):
    caps = si_captions()
    body = body_text()
    items = ["S1 Table", "S2 Table", "S3 Table", "S1 Text", "S2 Text", "S3 Text", "S1 Fig", "S2 Fig", "S3 Fig", "S4 Fig"]
    assert sorted(items) == sorted(caps), sorted(caps)
    files = {it: it.replace(" ", "_") + (".xlsx" if "Table" in it else ".docx" if "Text" in it else ".tif") for it in items}
    cites = {it: len(re.findall(r"\b%s\b" % re.escape(it), body)) for it in items}
    rows = []
    for it in items:
        if "Table" in it:
            rows.append([it, files[it], "index sheet", "Contents of this workbook.", caps[it][0], cites[it]])
            for sh in books[it]:
                rows.append([it, files[it], "%s sheet" % sh["name"], sh["title"][:200], caps[it][0], cites[it]])
        elif "Text" in it:
            src = (REPO / "docs" / "manuscript" / "si" / ("%s_source.txt" % files[it][:-5])).read_text()
            for m in re.finditer(r"^## (Section \d+)\. (.+)$", src, re.M):
                rows.append([it, files[it], m.group(1), m.group(2).strip(), caps[it][0], cites[it]])
        else:
            rows.append([it, files[it], "whole file", caps[it][0], caps[it][0], cites[it]])
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["item", "file", "locator", "title", "caption_title", "body_citations"])
    w.writerows(rows)
    index_csv = buf.getvalue()
    sizes = {}
    for it in items:
        p = SI / files[it]
        sizes[it] = round(p.stat().st_size / 1024)
    words = {10: "Ten"}
    lines = ["COMAVI supporting information", "=============================", "",
             "%s supporting-information items. Upload each file below to PLOS Computational Biology as a separate"
             % words[len(items)],
             "Supporting Information item, using the item name in the first column. The captions are listed at the end of",
             "the manuscript file, as PLOS requires; they are not repeated here and must not be uploaded separately.", "",
             "Supporting information is published exactly as provided and is not copyedited, so each file is",
             "self-describing: every workbook sheet opens on its own caption row and every document on a contents table.", "",
             "  item       file             what it contains",
             "  " + "-" * 94]
    for it in items:
        lines.append("  %-10s %-16s %s (%d KB)" % (it, files[it], caps[it][0], sizes[it]))
    lines += ["", "Items are numbered in the order the manuscript first cites them. The four figures are 600 dpi LZW TIFFs",
              "written by figures/submission/make_submission_figures.py and pinned in figures/submission/SHA256SUMS.txt, so",
              "any of them can be regenerated and checked. The three workbooks are written from the committed records by",
              "scripts/build_supplement.py, and the three texts from their sources by scripts/build_si_texts.py.",
              "SI_INDEX.csv lists every sheet and section with its title; it is a repository record and is not uploaded.", ""]
    readme = "\n".join(lines)
    sums = []
    for it in items:
        blob = (SI / files[it]).read_bytes()
        sums.append("%s  %s" % (hashlib.sha256(blob).hexdigest(), files[it]))
    sums.append("%s  SI_INDEX.csv" % hashlib.sha256(index_csv.encode()).hexdigest())
    return {"SI_INDEX.csv": index_csv, "README.txt": readme, "SHA256SUMS.txt": "\n".join(sums) + "\n"}


def workbook_content(data):
    """Every sheet as rows of cell values: what a reader sees, independent of
    the library version that serialized it."""
    wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True)
    out = [(name, [tuple(r) for r in wb[name].iter_rows(values_only=True)]) for name in wb.sheetnames]
    wb.close()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare; write nothing")
    args = ap.parse_args()
    SI = si_dir()
    books = {"S1 Table": s1_sheets(), "S2 Table": s2_sheets(), "S3 Table": s3_sheets()}
    caps = si_captions()
    for it in ("S1 Table", "S2 Table", "S3 Table"):
        assert caps[it][0] == TITLES[it[:2]], (it, caps[it][0])
    tables = {SI / ("%s.xlsx" % it.replace(" ", "_")): workbook_bytes(it, books[it])
              for it in ("S1 Table", "S2 Table", "S3 Table")}
    if args.check:
        # Workbooks are compared by content, since another openpyxl version
        # serializes the same cells to different bytes. The manifest is then
        # recomputed from the committed files, so SHA256SUMS.txt is checked
        # against what is actually shipped.
        bad = []
        for p, data in tables.items():
            if not p.exists() or workbook_content(p.read_bytes()) != workbook_content(data):
                bad.append("%s content differs from a fresh build" % p.relative_to(REPO))
        if not bad:
            for k, v in manifest(SI, books).items():
                p = SI / k
                if not p.exists() or p.read_text() != v:
                    bad.append("%s differs from the manifest of the committed files" % p.relative_to(REPO))
        for msg in bad:
            print("FAIL: %s" % msg)
        if bad:
            print("      run python scripts/build_supplement.py (after build_si_texts.py if a text changed)")
            return 1
        print("PASS: S1-S3 Table reproduce from the records (%d sheets), and SI_INDEX.csv, README.txt and "
              "SHA256SUMS.txt match the committed package" % sum(len(books[it]) for it in ("S1 Table", "S2 Table", "S3 Table")))
        return 0
    for p, data in tables.items():
        if not p.exists() or workbook_content(p.read_bytes()) != workbook_content(data):
            p.write_bytes(data)
            print("wrote %s (%d bytes)" % (p.relative_to(REPO), len(data)))
    for k, v in manifest(SI, books).items():
        (SI / k).write_text(v)
        print("wrote %s (%d bytes)" % ((SI / k).relative_to(REPO), len(v.encode())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
