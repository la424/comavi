#!/usr/bin/env python3
"""Recompute every derived cell of the supplementary workbooks.

WHY THIS EXISTS
---------------
Nothing checked the contents of the supplementary workbooks. The final v9
audit demonstrated this twice over. First, by finding a real defect the hard
way: the S3 Table tier_distribution sheet still read "Tier 2  6/17 = 0.353"
and asserted "Tiers 3-4  0/17 = 0.000" long after the ground truth moved, while
the S2 Text prose next to it already said 4/17. Second, by negative control:
editing a cell of S3_Table.xlsx to read "11/13 = 0.999" and re-running the
package, claims and denominator gates produced three passes. A checksum
manifest is not a substitute: the stale sheet was already stale when the
package was built, so its digest matched perfectly.

v11 rebuilt the workbooks from the records (scripts/build_supplement.py), and
this gate was rewritten with them. It does not import the builder. Each sheet
is recomputed here from the most primitive record that holds it -- the
canonical table, the evidence ledger, the per-variant prioritization table,
the calibration points and the SKEMPI tables -- and compared cell by cell;
values that exist only in a statistical record (bootstrap intervals,
permutation nulls) are compared with that record. A shown number passes when
it is within half a unit of its last displayed digit of the exact value.

Every run also plants a wrong cell in each checked sheet, in memory, and
requires the sheet's check to report it. A sheet whose check cannot see a
planted error fails the gate, so a check that silently stopped comparing
anything is caught.

Out of scope, and listed on every run: the index sheets (captions) and the
model_scope sheet, which is curated prose with no derived value.

Usage:
    PYTHONPATH=.:scripts python scripts/verify_si_derived_sheets.py
"""
from __future__ import annotations

import copy
import csv
import glob
import json
import pathlib
import re
import sys
from fractions import Fraction
from math import comb

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import fisher_exact, pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"
sys.path[:0] = [str(REPO / "scripts"), str(REPO)]
import benchmark_population as bp  # noqa: E402
import build_manuscript as bm  # noqa: E402

TAGS = ("t10", "t15", "t20", "t25", "tSAP")
GM = {"consistent": Fraction(1), "partial": Fraction(1, 2), "inconsistent": Fraction(0)}
CODE = {"consistent": "C", "partial": "P", "inconsistent": "X"}
OP_ROW = {"t10": "1.0 kcal/mol", "t15": "1.5 kcal/mol", "t20": "2.0 kcal/mol", "t25": "2.5 kcal/mol (reference)",
          "tSAP": "Prediction-interval limit (2.9 / 2.9 / 3.5)"}
SWEEP_ROW = {"t10": "1.0", "t15": "1.5", "t20": "2.0", "t25": "2.5 (reference)", "tSAP": "2.9/2.9/3.5",
             "t30": "3.0 (post hoc)", "t35": "3.5 (post hoc)", "t40": "4.0 (post hoc)"}
EXPECTED_SHEETS = {
    "S1_Table.xlsx": ["index", "evidence_ledger", "whole_variant_negatives"],
    "S2_Table.xlsx": ["index", "benchmark_variants", "evidence_composition", "model_scope", "excluded_variants",
                      "context_comparators", "fold_binding_disagreements"],
    "S3_Table.xlsx": ["index", "threshold_sweep", "operating_point_sweep", "missed_mechanisms", "disagreement_classes",
                      "tier_distribution", "evidence_states", "detection_metrics", "pathogenicity_ranking",
                      "robustness", "skempi_validation"],
}
EVIDENCE_LABEL = {"E1_quantitative_energetic": "E1 quantitative energetic",
                  "E2_quantitative_functional": "E2 quantitative functional",
                  "E3_qualitative_experimental": "E3 qualitative experimental",
                  "E4_population_frequency": "E4 population frequency",
                  "E5_inferred_no_axis_assay": "E5 inferred (no axis assay)"}
OUT_OF_SCOPE = {"index": "captions only", "model_scope": "curated prose; no derived value"}
NUM = re.compile(r"(?<![A-Za-z0-9])-?\d+(?:\.\d+)?")


def J(name):
    return json.loads((RO / name).read_text())


# ------------------------------------------------------------------ sheets
class Book:
    """Sheets as header + rows of strings (caption on row 1, header on row 3)."""

    def __init__(self, sheets):
        self.sheets = sheets

    @classmethod
    def read(cls, path):
        wb = openpyxl.load_workbook(path, read_only=True)
        out = {}
        for name in wb.sheetnames:
            rows = [["" if c is None else str(c) for c in r] for r in wb[name].iter_rows(values_only=True)]
            if name == "index":
                out[name] = {"header": rows[0], "rows": rows[1:]}
            else:
                out[name] = {"header": rows[2], "rows": [r for r in rows[3:] if any(x for x in r)]}
        wb.close()
        return cls(out)


class Ctx:
    """Collects comparisons and failures for one sheet."""

    def __init__(self, book, sheet, keycols=1):
        self.sh = book.sheets[sheet]
        self.name, self.k = sheet, keycols
        self.hdr = self.sh["header"]
        self.fail, self.cells, self.seen = [], [], set()
        self.index = {}
        for i, r in enumerate(self.sh["rows"]):
            key = tuple(r[:keycols]) if keycols > 1 else r[0]
            if key in self.index:
                self.fail.append("%s: duplicate row %s" % (sheet, key))
            self.index[key] = i

    def _cell(self, key, col):
        if key not in self.index:
            self.fail.append("%s: row %s missing" % (self.name, key))
            return None, None
        if col not in self.hdr:
            self.fail.append("%s: column %r missing" % (self.name, col))
            return None, None
        i, j = self.index[key], self.hdr.index(col)
        self.seen.add(key)
        return (i, j), self.sh["rows"][i][j]

    def num(self, key, col, *expected):
        """Numbers shown in the cell, in order, against exact expected values."""
        pos, shown = self._cell(key, col)
        if pos is None:
            return
        self.cells.append(pos)
        got = NUM.findall(shown.replace("\u2212", "-"))
        if len(got) != len(expected):
            self.fail.append("%s %s / %s: shows %r, expected %d numbers" % (self.name, key, col, shown, len(expected)))
            return
        for g, e in zip(got, expected):
            if e is None:
                continue
            if "." in g:
                d = len(g.split(".")[1])
                if abs(Fraction(g) - Fraction(e)) > Fraction(1, 2 * 10 ** d) + Fraction(1, 10 ** 9):
                    self.fail.append("%s %s / %s: shows %s, record %s" % (self.name, key, col, shown, float(e)))
                    return
            elif Fraction(int(g)) != Fraction(e):
                self.fail.append("%s %s / %s: shows %s, record %s" % (self.name, key, col, shown, float(e)))
                return

    def text(self, key, col, expected):
        pos, shown = self._cell(key, col)
        if pos is None:
            return
        if shown != expected:
            self.fail.append("%s %s / %s: shows %r, record %r" % (self.name, key, col, shown[:80], str(expected)[:80]))

    def rows_exactly(self, keys):
        extra = sorted(set(self.index) - set(keys), key=str)
        missing = sorted(set(keys) - set(self.index), key=str)
        if extra:
            self.fail.append("%s: rows not in the record: %s" % (self.name, extra[:6]))
        if missing:
            self.fail.append("%s: record rows missing from the sheet: %s" % (self.name, missing[:6]))


def F(x):
    """The exact rational behind a stored float (statistics of half-integer grades)."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, (int, np.integer)):
        return Fraction(int(x))
    f = Fraction(repr(float(x)))
    g = f.limit_denominator(10 ** 6)
    return g if abs(float(g) - float(x)) < 1e-12 else f


def none_if_dash(s):
    return None if s in ("", "-") else float(s)


# ------------------------------------------------------------------ records
CAN = pd.read_csv(RO / "scored_61var_canonical.csv", low_memory=False)
B = bp.benchmark(CAN).reset_index(drop=True)
MECH = B.expected_mech_class != "structurally_silent"
LED = pd.read_csv(RO / "COMAVI_evidence_ledger.csv", keep_default_na=False, dtype=str)
BK = set(zip(B.system, B.variant))
LB = LED[[k in BK for k in zip(LED.system, LED.variant)]]
CT = {(r["Gene"].lower(), r["Variant"]): r for r in
      csv.DictReader((RO / "cohort_tables" / "COMAVI_TableS1_variants.csv").open())}
PV = pd.read_csv(RO / "isds_v1" / "ISDS_v1_per_variant.csv")
NUMBERS = bm.number_citations(bm.parse_source(bm.SOURCE.read_text()))
KEYS = {r["evidence_citation"]: r["reference_keys"].split(";")
        for r in csv.DictReader((REPO / "data" / "analysis_inputs" / "si" / "ledger_citation_keys.csv").open())}


def grades(tag, mask=None):
    g = B["mech_consistency_%s" % tag].map(GM)
    return list(g if mask is None else g[mask])


def refnums(keys):
    return sorted(NUMBERS[k] for k in keys)


# ================================================================ S1 Table
def check_evidence_ledger(book):
    c = Ctx(book, "evidence_ledger", keycols=3)
    keys = []
    for r in LB.itertuples():
        k = (r.system, r.variant, r.axis)
        keys.append(k)
        for col in ("expected_token", "evidence_type", "evidence_directness", "evidence_basis", "evidence_citation"):
            c.text(k, col, getattr(r, col))
        c.text(k, "complex_amount_assay", r.complex_amount_assay or "-")
        c.num(k, "reference", *refnums(KEYS[r.evidence_citation]))
    c.rows_exactly(keys)
    return c


def check_whole_variant_negatives(book):
    c = Ctx(book, "whole_variant_negatives", keycols=3)
    want = sorted(BK - set(zip(LB.system, LB.variant)))
    rows = list(csv.DictReader((REPO / "data" / "analysis_inputs" / "si" / "whole_variant_negatives.csv").open()))
    gene = dict(zip(zip(B.system, B.variant), B.gene.astype(str).str.upper()))
    if sorted((d["system"], d["variant"]) for d in rows) != want:
        c.fail.append("whole_variant_negatives: curated input does not match the benchmark variants without a ledger row")
    keys = []
    for d in rows:
        k = (d["system"], gene[(d["system"], d["variant"])], d["variant"])
        keys.append(k)
        c.text(k, "evidence_basis", d["evidence_basis"])
        c.num(k, "reference", *refnums(d["reference_keys"].split(";")))
        row = B[(B.system == d["system"]) & (B.variant == d["variant"])].iloc[0]
        if row.expected_mech_class != "structurally_silent":
            c.fail.append("whole_variant_negatives: %s %s is not curated non-disruptive" % k[1:])
    c.rows_exactly(keys)
    return c


# ================================================================ S2 Table
def check_benchmark_variants(book):
    c = Ctx(book, "benchmark_variants", keycols=3)
    keys = []
    for _, r in B.iterrows():
        ct = CT[(str(r.gene).lower(), r.variant)]
        k = (ct["System"], str(r.gene).upper(), r.variant)
        keys.append(k)
        c.text(k, "Grade at 2.5 kcal/mol", r.mech_consistency_t25)
        c.text(k, "COMAVI call at 2.5 kcal/mol", r.mech_t25)
        codes = "/".join(CODE[r["mech_consistency_%s" % t]] for t in TAGS)
        c.text(k, "Grades at the five operating points", codes)
        c.text(k, "Same grade at every operating point", "yes" if len(set(codes.split("/"))) == 1 else "no")
        # the BRCT system is a single chain, which the tier does not score
        c.text(k, "Structural tier", "-" if pd.isna(r.comavi_tier) else r.comavi_tier)
        for col, src in (("Predicted \u0394\u0394G, isolated subunit", "Predicted DDG, isolated subunit (kcal/mol)"),
                         ("Predicted \u0394\u0394G, assembled complex", "Predicted DDG, in complex (kcal/mol)"),
                         ("Predicted \u0394\u0394G, binding", "Predicted DDG, binding (kcal/mol)"),
                         ("Measured \u0394\u0394G", "Measured DDG (kcal/mol)")):
            v = none_if_dash(ct[src])
            if v is None:
                c.text(k, col, "-")
            else:
                c.num(k, col, v)
        am = r.get("AM pathogenicity")
        if pd.isna(am):
            c.text(k, "AlphaMissense pathogenicity", "n/a")
        else:
            c.num(k, "AlphaMissense pathogenicity", am)
    c.rows_exactly(keys)
    return c


def check_evidence_composition(book):
    c = Ctx(book, "evidence_composition")
    for et in sorted(LB.evidence_type.unique()):
        sub = LB[LB.evidence_type == et]
        key = EVIDENCE_LABEL[et]
        tok, d = sub.expected_token.value_counts(), sub.evidence_directness.value_counts()
        for col, v in (("Expectations", len(sub)), ("Neutral", tok.get("neutral", 0)), ("Destabilizing", tok.get("destab", 0)),
                       ("Stabilizing", tok.get("stab", 0)), ("Direct", d.get("direct", 0)), ("Coupled", d.get("coupled", 0)),
                       ("Inferred", d.get("inferred", 0))):
            c.num(key, col, int(v))
    c.num("All classes", "Expectations", len(LB))
    return c


def check_excluded_variants(book):
    c = Ctx(book, "excluded_variants", keycols=3)
    ex = bp.excluded(CAN)
    keys = []
    for _, r in ex.iterrows():
        ct = CT[(str(r.gene).lower(), r.variant)]
        k = (ct["System"], str(r.gene).upper(), r.variant)
        keys.append(k)
        c.text(k, "Reason", r.exclusion_reason)
        c.text(k, "COMAVI call at 2.5 kcal/mol", r.mech_t25)
        for col, src in (("Predicted \u0394\u0394G, isolated subunit", "Predicted DDG, isolated subunit (kcal/mol)"),
                         ("Predicted \u0394\u0394G, assembled complex", "Predicted DDG, in complex (kcal/mol)"),
                         ("Predicted \u0394\u0394G, binding", "Predicted DDG, binding (kcal/mol)")):
            v = none_if_dash(ct[src])
            if v is None:
                c.text(k, col, "-")
            else:
                c.num(k, col, v)
    c.rows_exactly(keys)
    if len(B) + len(keys) != len(CAN):
        c.fail.append("excluded_variants: benchmark + excluded != scored cohort")
    return c


def perm_p(system, truth, strong):
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


def check_context_comparators(book):
    c = Ctx(book, "context_comparators")
    y = PV.structural_ground_truth.astype(bool)
    for name, strong in (("Interface status alone (any confidence-passing interface contact)",
                          PV.interface_only_positive.astype(bool)),
                         ("Tier 1-2 without the interface bonus", PV.no_interface_score >= 3.0),
                         ("Tier 1-2 (full tier)", PV.tier_number.isin([1, 2]))):
        tp, fn = int((strong & y).sum()), int((~strong & y).sum())
        fp, tn = int((strong & ~y).sum()), int((~strong & ~y).sum())
        c.num(name, "Mechanisms recovered / missed", tp, fn)
        c.num(name, "Non-disruptive called / rejected", fp, tn)
        c.num(name, "Sensitivity", Fraction(tp, tp + fn))
        c.num(name, "Specificity", Fraction(tn, tn + fp))
        c.num(name, "Positive predictive value", Fraction(tp, tp + fp))
        c.num(name, "Balanced accuracy", (Fraction(tp, tp + fn) + Fraction(tn, tn + fp)) / 2)
        c.num(name, "Fisher exact p", F(fisher_exact([[tp, fn], [fp, tn]])[1]))
        c.num(name, "Within-system permutation p", perm_p(PV.system, y, strong))
    return c


def check_fold_binding_disagreements(book):
    c = Ctx(book, "fold_binding_disagreements")
    cal = dict(zip(*[pd.read_csv(RO / "COMAVI_delta_calibration_points.csv")[k] for k in ("variant", "measured_kcal")]))
    col = {"Isolated-subunit stability": "Predicted DDG, isolated subunit (kcal/mol)",
           "Binding": "Predicted DDG, binding (kcal/mol)"}
    for r in book.sheets["fold_binding_disagreements"]["rows"]:
        name, axis = r[0], r[1]
        gene, var = name.split()
        ct = CT[(gene.lower(), var)]
        if (CAN[(CAN.gene.astype(str).str.upper() == gene) & (CAN.variant == var)].system.iloc[0], var) not in BK:
            c.fail.append("fold_binding_disagreements: %s is not a benchmark variant" % name)
        if axis in col:
            c.num(name, "Predicted \u0394\u0394G", float(ct[col[axis]]))
            shown = r[2]
            if NUM.search(shown):
                c.num(name, "Measured or curated expectation", float(cal[var]))
        else:
            c.num(name, "Predicted \u0394\u0394G", float(ct["Predicted DDG, isolated subunit (kcal/mol)"]),
                  float(ct["Predicted DDG, in complex (kcal/mol)"]))
    return c


# ================================================================ S3 Table
def check_threshold_sweep(book):
    c = Ctx(book, "threshold_sweep")
    sw = {r["threshold_tag"]: r for r in J("COMAVI_threshold_sweep.json")["rows"]}
    ext = {r["threshold_tag"]: r for r in J("COMAVI_threshold_sweep_extended.json")["rows"]}
    cols = [("Detection", "detection"), ("Attribution", "attribution"), ("Correct rejection", "correct_rejection"),
            ("Axis: isolated subunit", "axis_monomer"), ("Axis: assembled complex", "axis_fold"),
            ("Axis: binding", "axis_binding"), ("Axis: either assembled axis (amount of complex)", "axis_assembled_either"),
            ("Detection, isolated subunit only", None), ("Correct rejection, isolated subunit only", None)]
    for t, label in SWEEP_ROW.items():
        rec = sw.get(t) or ext[t]
        for col, stem in cols:
            if stem is None:
                k = rec["detection_isolated_k" if col.startswith("Detection") else "correct_rejection_isolated_k"]
                n = rec["detection_n" if col.startswith("Detection") else "correct_rejection_n"]
            else:
                k, n = rec[stem + "_k"], rec[stem + "_n"]
            c.num(label, col, Fraction(int(k), int(n)))
        if t in TAGS:
            # recomputed from the canonical table, not the sweep record
            c.num(label, "Whole-variant score", sum(grades(t)) / len(B))
            c.num(label, "Structural agreement", Fraction(int(B["structural_agreement_n_%s" % t].sum()),
                                                          int(B["structural_agreement_d_%s" % t].sum())))
        else:
            c.text(label, "Whole-variant score", "-")
            c.num(label, "Structural agreement", Fraction(int(rec["structural_agreement_k"]), int(rec["structural_agreement_n"])))
    c.rows_exactly(list(SWEEP_ROW.values()))
    return c


def check_operating_point_sweep(book):
    c = Ctx(book, "operating_point_sweep")
    opd = J("COMAVI_operating_point_diagnostics.json")["paired_bootstrap"]
    rec = {d["threshold"]: d for d in J("COMAVI_measured_effect_tables.json")["measured_effects_recovered"]["by_threshold"]}
    for t, label in OP_ROW.items():
        am = sum(grades(t, MECH)) / int(MECH.sum())
        asl = sum(grades(t, ~MECH)) / int((~MECH).sum())
        c.num(label, "Mechanism arm (n = 29)", am, *map(F, opd["recovery"][t]["ci95"]))
        c.num(label, "Non-disruptive arm (n = 27)", asl, *map(F, opd["correct_rejection"][t]["ci95"]))
        c.num(label, "Pooled score (n = 56)", sum(grades(t)), len(B), sum(grades(t)) / len(B))
        c.num(label, "Balanced score", (am + asl) / 2)
        m = rec[{"t10": "1.0", "t15": "1.5", "t20": "2.0", "t25": "2.5", "tSAP": "tSAP"}[t]]
        c.num(label, "Measured destabilizations recovered (n = 42)", m["recovered"], m["n"])
    c.rows_exactly(list(OP_ROW.values()))
    return c


def check_missed_mechanisms(book):
    c = Ctx(book, "missed_mechanisms", keycols=2)
    msd = J("COMAVI_detection_vs_attribution.json")["missed_at_reference"]
    keys = []
    for d in msd["variants"]:
        k = (d["gene"], d["variant"])
        keys.append(k)
        r = CAN[(CAN.gene.astype(str).str.upper() == d["gene"]) & (CAN.variant == d["variant"])].iloc[0]
        c.text(k, "COMAVI call at 2.5 kcal/mol", r.mech_t25)
        for col, f in (("Predicted \u0394\u0394G, isolated subunit", "predicted_monomer"),
                       ("Predicted \u0394\u0394G, assembled complex", "predicted_complex_fold"),
                       ("Predicted \u0394\u0394G, binding", "predicted_binding")):
            if d[f] is None:
                c.text(k, col, "-")
            else:
                c.num(k, col, d[f])
        if d.get("measured_kcal_mol") is None:
            c.text(k, "Measured \u0394\u0394G", "-")
        else:
            c.num(k, "Measured \u0394\u0394G", d["measured_kcal_mol"])
            c.text(k, "Measured effect", "below the reference" if d["measured_kcal_mol"] < 2.5 else "above the reference")
    c.rows_exactly(keys)
    # the set itself: structural mechanisms of the benchmark whose t25 grade has no detected axis
    struct = set(zip(B[MECH].gene.astype(str).str.upper(), B[MECH].variant))
    if not set(keys) <= struct:
        c.fail.append("missed_mechanisms: lists variants that are not benchmark structural mechanisms")
    return c


def check_disagreement_classes(book):
    c = Ctx(book, "disagreement_classes")
    cal = pd.read_csv(RO / "COMAVI_delta_calibration_points.csv")
    st = cal.straddles.fillna(False).astype(bool)
    near = cal.both_near.fillna(False).astype(bool)
    err = (cal.foldx_ddg - cal.measured_kcal).abs()
    g73 = cal.g73.fillna(False).astype(bool)
    for name, mask in (("Agreement", ~st), ("Borderline disagreement", st & near), ("Substantive disagreement", st & ~near)):
        c.num(name, "n", int(mask.sum()))
        c.num(name, "Median absolute error (kcal/mol)", F(float(err[mask].median())))
        c.num(name, "n without Glu73", int((mask & ~g73).sum()))
        c.num(name, "Median without Glu73 (kcal/mol)", F(float(err[mask & ~g73].median())))
    return c


def check_tier_distribution(book):
    c = Ctx(book, "tier_distribution")
    y = PV.structural_ground_truth.astype(bool)
    for label, mask in (("Tier 1", PV.tier_number == 1), ("Tier 2", PV.tier_number == 2), ("Tier 3", PV.tier_number == 3),
                        ("Tier 4", PV.tier_number == 4), ("Tiers 3-4", PV.tier_number.isin([3, 4]))):
        k, n = int((mask & y).sum()), int(mask.sum())
        c.num(label, "Structural mechanism", k)
        c.num(label, "Non-disruptive", n - k)
        c.num(label, "Structural fraction", k, n, Fraction(k, n))
    return c


def check_evidence_states(book):
    c = Ctx(book, "evidence_states")
    for row in J("isds_v1/ISDS_v1_evidence_state_tables.json")["four_evidence_states"]:
        for s in ("Convergent", "Energy only", "Context only", "Neither"):
            c.num(OP_ROW[row["threshold"]], s, row[s]["structural"], row[s]["no_lesion"])
    return c


def check_detection_metrics(book):
    c = Ctx(book, "detection_metrics")
    for row in J("isds_v1/ISDS_v1_evidence_state_tables.json")["detection_metrics"]:
        for key, lab in (("any_energetic_firing", "Any energy call"), ("convergent_energy_and_tier", "Convergent")):
            q = row[key]
            c.num(OP_ROW[row["threshold"]], lab + ": counts", q["tp"], q["fn"], q["fp"], q["tn"])
            c.num(OP_ROW[row["threshold"]], lab + ": rates", Fraction(q["tp"], q["tp"] + q["fn"]),
                  Fraction(q["tn"], q["tn"] + q["fp"]), Fraction(q["tp"], q["tp"] + q["fp"]))
    return c


def check_pathogenicity_ranking(book):
    c = Ctx(book, "pathogenicity_ranking")
    cs = pd.read_csv(RO / "isds_v1" / "ISDS_v1_alphamissense_common_set.csv")
    y = cs.clinical_y.astype(int)
    n1, n0 = int(y.sum()), int((1 - y).sum())
    for label, col in (("AlphaMissense", "AM pathogenicity"), ("Priority score", "isds_v1"),
                       ("Energy component of the priority score", "isds_energy_component"),
                       ("Context component of the priority score", "isds_context_component")):
        auc = roc_auc_score(y, cs[col])
        c.num(label, "Variants (pathogenic / benign)", len(cs), n1, n0)
        c.num(label, "ROC AUC for pathogenicity", Fraction(int(round(auc * n1 * n0 * 2)), 2 * n1 * n0))
    return c


def check_robustness(book):
    c = Ctx(book, "robustness")
    z = np.load(RO / "stress_tests" / "comavi_stress_draws.npz")
    vb = J("stress_tests/COMAVI_stress_verified_statistics.json")
    tot = sum(grades("t25"))
    sk, sn = int(B.structural_agreement_n_t25.sum()), int(B.structural_agreement_d_t25.sum())
    mc, sa = tot / len(B), Fraction(sk, sn)
    c.num("Observed at 2.5 kcal/mol", "Mechanism-pattern score", tot, len(B), mc)
    c.num("Observed at 2.5 kcal/mol", "Structural agreement", sk, sn, sa)
    lab = "Label-permutation null, %s permutations: mean (SD); p" % format(len(z["mc_null"]), ",")
    for col, null, obs in (("Mechanism-pattern score", z["mc_null"], mc), ("Structural agreement", z["sa_null"], sa)):
        p = Fraction(int((null >= float(obs) - 1e-12).sum()) + 1, len(null) + 1)
        c.num(lab, col, F(float(null.mean())), F(float(null.std())), p)
    g = B.assign(g=grades("t25"), k=B.structural_agreement_n_t25, n=B.structural_agreement_d_t25).groupby("system")
    mcs = [(tot - Fraction(gg.g.sum())) / (len(B) - len(gg)) for _, gg in g]
    sas = [Fraction(sk - int(gg.k.sum()), sn - int(gg.n.sum())) for _, gg in g]
    lo = "Leave one system out (%d systems): range" % g.ngroups
    c.num(lo, "Mechanism-pattern score", min(mcs), max(mcs))
    c.num(lo, "Structural agreement", min(sas), max(sas))
    nz = "Replicate noise, %d draws: mean (SD)" % len(z["noise_mc"])
    c.num(nz, "Mechanism-pattern score", F(float(z["noise_mc"].mean())), F(float(z["noise_mc"].std())))
    c.num(nz, "Structural agreement", F(float(z["noise_sa"].mean())), F(float(z["noise_sa"].std())))
    for lab2, key, dk in (("Variant bootstrap", "variant_bootstrap_ci95", "variant_bootstrap"),
                          ("System-cluster bootstrap", "system_cluster_ci95", "system_cluster_bootstrap")):
        row = "%s, %s draws: 95%% interval" % (lab2, format(vb["draws"][dk], ","))
        c.num(row, "Mechanism-pattern score", *map(F, vb[key]["mc"]["full_precision"]))
        c.num(row, "Structural agreement", *map(F, vb[key]["sa"]["full_precision"]))
    return c


def check_skempi_validation(book):
    c = Ctx(book, "skempi_validation")
    bb = pd.read_csv(REPO / "supplement" / "skempi" / "bb_binding_validation.csv")
    jt = pd.read_csv(REPO / "supplement" / "skempi" / "jt_binding_validation.csv")
    g73 = bb.foldx_code.str.match(r"E[A-Z]?73")
    for label, df in (("Barnase-barstar, all", bb), ("Barnase-barstar, without Glu73", bb[~g73]),
                      ("TEM1-BLIP, BLIP side", jt), ("Both complexes", pd.concat([bb, jt])),
                      ("Both complexes, without Glu73", pd.concat([bb[~g73], jt]))):
        r, pr = pearsonr(df.ddg_meas, df.ddg_pred_bind)
        rho, prho = spearmanr(df.ddg_meas, df.ddg_pred_bind)
        d = df[df.ddg_meas > 0]
        c.num(label, "Substitutions", len(df))
        c.num(label, "Pearson r", F(float(r)), F(float(pr)))
        c.num(label, "Spearman \u03c1", F(float(rho)), F(float(prho)))
        c.num(label, "Sign correct on measured destabilizations", int((d.ddg_pred_bind > 0).sum()), len(d))
    return c


CHECKS = {
    "S1_Table.xlsx": {"evidence_ledger": check_evidence_ledger, "whole_variant_negatives": check_whole_variant_negatives},
    "S2_Table.xlsx": {"benchmark_variants": check_benchmark_variants, "evidence_composition": check_evidence_composition,
                      "excluded_variants": check_excluded_variants, "context_comparators": check_context_comparators,
                      "fold_binding_disagreements": check_fold_binding_disagreements},
    "S3_Table.xlsx": {"threshold_sweep": check_threshold_sweep, "operating_point_sweep": check_operating_point_sweep,
                      "missed_mechanisms": check_missed_mechanisms, "disagreement_classes": check_disagreement_classes,
                      "tier_distribution": check_tier_distribution, "evidence_states": check_evidence_states,
                      "detection_metrics": check_detection_metrics, "pathogenicity_ranking": check_pathogenicity_ranking,
                      "robustness": check_robustness, "skempi_validation": check_skempi_validation},
}


def plant(book, sheet, pos):
    """A copy of the book with the first number of one checked cell moved by one unit of 10."""
    b = copy.deepcopy(book)
    i, j = pos
    cell = b.sheets[sheet]["rows"][i][j]
    m = NUM.search(cell)
    v = m.group(0)
    wrong = str(int(v) + 10) if "." not in v else "%.*f" % (len(v.split(".")[1]), float(v) + 10)
    b.sheets[sheet]["rows"][i][j] = cell[:m.start()] + wrong + cell[m.end():]
    return b


def si_dir():
    hits = sorted(glob.glob(str(REPO / "submission" / "COMAVI_*_SI")))
    if len(hits) != 1:
        raise SystemExit("FAIL: expected exactly one submission/COMAVI_*_SI directory, found %s" % hits)
    return pathlib.Path(hits[0])


def main() -> int:
    si = si_dir()
    failures, n_cells, controls = [], 0, []
    for fname, checks in CHECKS.items():
        book = Book.read(si / fname)
        if list(book.sheets) != EXPECTED_SHEETS[fname]:
            failures.append("%s sheets %s, expected %s" % (fname, list(book.sheets), EXPECTED_SHEETS[fname]))
        for sheet in EXPECTED_SHEETS[fname]:
            if sheet not in checks and sheet not in OUT_OF_SCOPE:
                failures.append("%s/%s has no check and is not declared out of scope" % (fname, sheet))
        for sheet, fn in checks.items():
            if sheet not in book.sheets:
                continue
            ctx = fn(book)
            failures += ctx.fail
            n_cells += len(ctx.cells)
            if not ctx.cells:
                failures.append("%s/%s: no numeric cell was compared" % (fname, sheet))
                continue
            # negative control: one planted wrong number must be reported
            planted = fn(plant(book, sheet, ctx.cells[0]))
            caught = len(planted.fail) > len(ctx.fail)
            controls.append((sheet, caught))
            if not caught:
                failures.append("%s/%s: a planted wrong cell was not caught" % (fname, sheet))
    if failures:
        print("FAIL: %d supplementary cells or structures disagree with their record" % len(failures))
        for f in failures[:60]:
            print("  " + f)
        return 1
    print("PASS: %d derived cells in %d sheets of S1-S3 Table recompute from their records; a planted wrong cell "
          "was caught in each of the %d sheets" % (n_cells, len(controls), sum(c for _, c in controls)))
    print("  out of scope: %s" % "; ".join("%s (%s)" % kv for kv in OUT_OF_SCOPE.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
