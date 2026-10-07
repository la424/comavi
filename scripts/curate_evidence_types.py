#!/usr/bin/env python
"""Per-axis evidence-type ledger for the COMAVI benchmark.

WHY THIS EXISTS
---------------
Every COMAVI metric is scored against a per-axis ground-truth token
(`expected_ddg_{monomer,fold_complex,binding}` in {destab, neutral, stab}).
Until now the paper could not report *what kind of evidence* backs each token,
so a reader could not tell a directly measured free energy from a token
inferred from a variant's position. That is a fair question and this file
answers it, at the resolution the metrics are actually computed on: the axis,
not the variant.

TWO ORTHOGONAL DESCRIPTORS
--------------------------
Conflating "is it quantitative?" with "does it speak to THIS axis?" is the
trap. A Ca2+-affinity K_d is a hard number that says nothing about an IQ-motif
interface; a co-IP "abolished" is qualitative but lands squarely on the binding
axis. So each axis gets both:

evidence_type      — the measurement's nature
  E1_quantitative_energetic ..... a free energy / affinity / thermal stability
                                  (dG_unfold, assembly ddG, K_D, Tm)
  E2_quantitative_functional .... a numeric non-energetic assay
                                  (co-IP % of WT, EC50, MAVE score, MMR %)
  E3_qualitative_experimental ... a categorical experimental result
                                  (abolished / preserved / natively folded)
  E4_population_frequency ....... benign token from allele frequency alone
  E5_inferred_no_axis_assay ..... no assay addresses this axis; the token is
                                  inferred from position, chemistry, or the
                                  absence of a reported effect

evidence_directness — what the measurement was performed on
  direct ...... the assay reads out this axis
  coupled ..... the assay reads out a process that REQUIRES this axis
                (MMR proficiency implies the MLH1-PMS2 interaction survived);
                the ledger's "relaxed" promotions live here
  off_axis .... quantitative, but on a different physical process than any
                structural axis (the ledger's reserved FUNCTIONAL tier)
  inferred .... no measurement; judgment from structure or population data

The pair is what makes the claim auditable. "E1/direct" is the strongest cell
and is what the paper's physical-validation sections rest on; "E5/inferred" is
the weakest and its prevalence is a real limitation the paper should state
rather than let a reviewer discover.

SCOPE
-----
Only axes carrying a COMMITTED token are curated. `unknown` and
`not_applicable` axes are excluded by construction: they are not scored, so
they have no ground truth to characterize. This is why the denominator here
(committed axes) differs from every other denominator in the paper -- see
docs/COMAVI_denominator_reconciliation.md.

Emits reference_outputs/COMAVI_evidence_ledger.csv (one row per committed
axis) and a summary JSON. Curation is asserted complete: every committed axis
must appear in ASSIGN or the script fails.
"""
from __future__ import annotations

import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
BRCT = REPO / "supplement" / "brct" / "brct_foldx_concordance.csv"
OUT_CSV = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"
OUT_JSON = REPO / "reference_outputs" / "COMAVI_evidence_ledger_summary.json"

AXES = {
    "monomer": "expected_ddg_monomer",
    "fold_complex": "expected_ddg_fold_complex",
    "binding": "expected_ddg_binding",
}
COMMITTED = {"destab", "neutral", "stab"}

E1, E2, E3, E4, E5 = (
    "E1_quantitative_energetic",
    "E2_quantitative_functional",
    "E3_qualitative_experimental",
    "E4_population_frequency",
    "E5_inferred_no_axis_assay",
)

# Citation shorthands, expanded once so a fix propagates.
C_STARITA = "Starita 2015 (Y2H interaction panel)"
C_BRZOVIC = ("Brzovic 1998 JBC 273:7795 (PMID 9520298); "
             "Roehm & Berg 1997 Biochemistry 36:10240")
C_CLARK = "Clark 2022 (mammalian two-hybrid)"
C_ZHAO = "Zhao & Vogt 2008 PNAS"
C_MANDELKER = "Mandelker 2009"
C_WU = "Wu / Jaiswal 2009"
C_KOSINSKI = "Kosinski 2010 (PMS2 co-IP, MMR activity)"
C_MAVE = "MAVE abundance/interaction scores; Kosinski 2010 co-IP"
C_OLLILA = ("Ollila 2008 Hum Mutat 29:1355 (PMID 18951462); "
            "Ollila 2006 Gastroenterology (PMID 17101317)")
C_JIA = "Jia 2021"
C_STEBBINS = "Stebbins 1999 Science"
C_OHH = "Ohh 2000"
C_HUNTER = ("Hunter 2015 Mol Cancer Res 13:1325 Table 2 "
            "(doi:10.1158/1541-7786.MCR-15-0203)")
C_HINRICHSEN = ("Hinrichsen 2013 Clin Cancer Res 19:2432 "
                "(doi:10.1158/1078-0432.CCR-12-3299)")
C_LINDSAY = ("Lindsay 2025 J Allergy Clin Immunol Fig 4D-E "
             "(doi:10.1016/j.jaci.2024.08.024)")
C_EATON = "Eaton & Hofrichter 1990 (HbS polymerization)"
C_KIGER = ("Kiger / Kwiatkowski Biochemistry 1998 "
           "(PMID 9521754 equilibrium assembly dG; PMID 9521753 kinetics)")
C_BONAVENTURA = "Bonaventura & Riggs 1968 JBC (Hb Kansas)"
C_ELLIOTT = ("Elliott 2000 JBC 275:22069; "
             "Takahashi / Kimura 2001 (PMID 11735257)")
C_KIMURA = "Kimura 2001 (PMID 11735257)"
C_FRAZIER = "Frazier 2008; gnomAD"
C_CROTTI = "Crotti 2013"
C_HWANG = "Hwang 2014 (PMC5270481)"
C_DEBOSSCHER = "De Bosscher 2004 (co-IP)"
C_LANAUZE = "Lanauze 2021 (co-IP)"
C_LEGOFF = "Le Goff 2011"
C_ROWLING = ("Rowling, Cook & Itzhaki 2010 JBC 285:20080 "
             "(PMID 20452977; PMC2888420), Table S1")
C_PECHTL = "Pechtl 2011 JBC 286:11082 (PMID 21270465), SPR"
C_TISCHER = ("Tischer / Moon-Tasson / Auton 2025 J Thromb Haemost "
             "23(4):1215-1228 (PMID 39756657); CD/DSC + SPR + shear-flow")
C_INSIGHT = "InSiGHT expert-panel reclassification; rs63750447 AF ~8.9% EAS"
C_GNOMAD = "gnomAD; ClinVar BA1"
C_RAEVAARA = "Raevaara 2005 (MMR activity)"
C_DROST = "Drost 2019 (functional assay)"

NO_STAB_ASSAY = "No subunit-stability assay; token inferred"
NO_ASSEMBLY_ASSAY = ("No assembly-context stability measurement; token "
                     "inferred from the intact monomer fold")

# (system, variant, axis) -> (type, directness, basis)
ASSIGN: dict[tuple[str, str, str], tuple[str, str, str]] = {}


def put(system, variant, axis, etype, directness, basis, citation):
    ASSIGN[(system, variant, axis)] = (etype, directness, basis, citation)


# ---------------------------------------------------------------- BRCA1-BARD1
put("brca1_bard1", "V11G", "binding", E2, "direct",
    "Moderate BARD1 binding loss by two-hybrid", C_CLARK)
for ax, basis in [
    ("monomer", "C61 is a Zn-coordinating ligand; RING fold lost on Zn release"),
    ("fold_complex", "RING heterodimer fold disrupted by Zn-coordination loss"),
]:
    put("brca1_bard1", "C61G", ax, E3, "direct", basis, C_BRZOVIC)
put("brca1_bard1", "C61G", "binding", E3, "direct",
    "BARD1 binding preserved despite fold loss", C_BRZOVIC)
# Values read from the source's own Supp_TABLE1 (5,160 variants; Y2H_score is
# the BARD1 two-hybrid interaction readout, 1.0 = wild type). Among the 1,472
# variants passing the authors' 800-read filter the median Y2H score is 0.963
# and the 5th percentile is 0.103, so both values below sit squarely in the
# wild-type part of the distribution.
#
# K45R  Y2H 1.019 (95% CI 0.983-1.056), filter Pass.
# E100Q Y2H 1.025 (95% CI 0.840-1.213), permutation p = 0.579, filter DROP.
#
# The E100Q measurement does not meet the authors' own read-depth filter and
# its interval is four times wider than K45R's. The point estimate is still
# wild-type-like with no evidence of a deficit, so the neutral token stands,
# but the basis records the filter status rather than presenting the two as
# equally supported -- a measurement the source itself declined to report is
# weaker evidence than one it reports.
for v in ("K45R", "E100Q"):
    put("brca1_bard1", v, "monomer", E5, "inferred",
        f"{NO_STAB_ASSAY} from WT-like two-hybrid function", C_STARITA)
    put("brca1_bard1", v, "fold_complex", E5, "inferred",
        NO_ASSEMBLY_ASSAY, C_STARITA)
put("brca1_bard1", "K45R", "binding", E3, "direct",
    "WT-like BARD1 interaction by Y2H: score 1.019 (95% CI 0.983-1.056), "
    "passes the source's 800-read filter", C_STARITA)
put("brca1_bard1", "E100Q", "binding", E3, "direct",
    "WT-like BARD1 interaction by Y2H: score 1.025 (95% CI 0.840-1.213), "
    "permutation p = 0.58, but the measurement is below the source's own "
    "800-read filter and is reported here with that caveat", C_STARITA)

# ----------------------------------------------------------------------- PI3K
for v in ("E542K", "E545K"):
    put("pi3k", v, "monomer", E5, "inferred",
        f"{NO_STAB_ASSAY}; helical-domain hotspot taken as fold-neutral", C_ZHAO)
    for ax in ("fold_complex", "binding"):
        put("pi3k", v, ax, E3, "direct",
            "Disrupts the nSH2 inhibitory contact, releasing autoinhibition",
            C_ZHAO)
put("pi3k", "H1047R", "monomer", E5, "inferred", NO_STAB_ASSAY, C_MANDELKER)
put("pi3k", "H1047R", "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY,
    C_MANDELKER)
put("pi3k", "H1047R", "binding", E3, "direct",
    "Kinase C-lobe gain of function with p85 binding structurally unaffected",
    C_MANDELKER)
put("pi3k", "N564D", "monomer", E5, "inferred", NO_STAB_ASSAY, C_WU)
put("pi3k", "N564D", "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY, C_WU)
put("pi3k", "N564D", "binding", E3, "direct",
    "Retains p110 binding, loses inhibition -- allosteric, not interface loss",
    C_WU)
for ax in AXES:
    put("pi3k", "M326I", ax, E4, "inferred",
        "rs3730089, gnomAD MAF > 10%: benign by frequency alone",
        "gnomAD (rs3730089)")

# ------------------------------------------------------------------ MLH1-PMS2
put("mlh1_pms2", "Q542L", "binding", E2, "direct",
    "PMS2 co-IP 54% of WT", C_KOSINSKI)
put("mlh1_pms2", "L749P", "fold_complex", E5, "inferred",
    "Proline substitution presumed fold-destabilizing; no stability assay",
    C_KOSINSKI)
put("mlh1_pms2", "L749P", "binding", E2, "direct",
    "PMS2 co-IP 30% of WT", C_KOSINSKI)
for ax in ("monomer", "fold_complex"):
    put("mlh1_pms2", "R755W", ax, E2, "direct",
        "MAVE abundance score -0.334: reduced cellular abundance", C_MAVE)
put("mlh1_pms2", "R755W", "binding", E2, "direct",
    "MAVE interaction score +0.617; PMS2 co-IP 82% -- binding preserved",
    C_MAVE)
put("mlh1_pms2", "V384D", "monomer", E5, "inferred", NO_STAB_ASSAY, C_INSIGHT)
put("mlh1_pms2", "V384D", "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY,
    C_INSIGHT)
put("mlh1_pms2", "V384D", "binding", E3, "coupled",
    "Reclassified non-pathogenic; MMR-proficient", C_INSIGHT)
for ax in AXES:
    put("mlh1_pms2", "G857A", ax, E4, "inferred",
        "gnomAD AF 28%; ClinVar BA1 stand-alone benign", C_GNOMAD)
# v8.0: a stability measurement for this exact variant does exist. Hinrichsen
# 2013 states that "the 4 MMR-proficient variants (K618A, E578G, V716M, and
# H718Y) moderately destabilize the MLH1 protein but are not causative for
# Lynch syndrome", with thermal stability "decreased more strongly in V716M,
# H718Y, and E578G" and pulse-chase half-life reduced to 64% on average.
# Structurally destabilizing and clinically neutral in the same variant.
put("mlh1_pms2", "H718Y", "monomer", E2, "direct",
    "Source states H718Y moderately destabilizes the MLH1 protein; thermal "
    "stability decreased more strongly in V716M, H718Y and E578G, and "
    "pulse-chase half-life fell to 64% of wild type. Clinically neutral",
    C_HINRICHSEN)
put("mlh1_pms2", "H718Y", "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY,
    C_RAEVAARA)
put("mlh1_pms2", "H718Y", "binding", E2, "coupled",
    "Proficient MMR implies the MLH1-PMS2 interaction survived", C_RAEVAARA)
put("mlh1_pms2", "K618E", "monomer", E5, "inferred", NO_STAB_ASSAY, C_DROST)
put("mlh1_pms2", "K618E", "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY,
    C_DROST)
put("mlh1_pms2", "K618E", "binding", E3, "coupled",
    "Non-damaging in functional assay", C_DROST)

# ------------------------------------------------------------------ MSH2-MSH6
# v8.0: G674R's binding commitment is retired, and the canonical now marks this
# axis `unknown`, so this entry is deliberately absent rather than corrected.
# The prior basis read "MSH6 binding preserved by co-IP" and cited Jia 2021,
# which measures MMR loss of function by 6-thioguanine selection and performs
# no co-immunoprecipitation. Ollila 2008 did run co-IP on an MSH2 ATPase-domain
# variant, but on G674A -- a different substitution at the same position; that
# paper never mentions G674R. No accessible source gives this variant's binding
# direction, so there is nothing to commit.
put("msh2_msh6", "A636P", "monomer", E3, "direct",
    "MSH2 expression/stability preserved; ATPase mismatch binding/release defect",
    C_OLLILA)
put("msh2_msh6", "A636P", "fold_complex", E3, "direct",
    "MSH2-MSH6 interaction intact; ATPase mismatch binding/release defect",
    C_OLLILA)
put("msh2_msh6", "A636P", "binding", E3, "direct",
    "MSH2-MSH6 interaction intact; ATPase mismatch binding/release defect",
    C_OLLILA)
put("msh2_msh6", "C697F", "binding", E3, "direct",
    "ATPase catalytic defect, protein stable, interaction intact", C_OLLILA)
put("msh2_msh6", "N127S", "binding", E3, "direct",
    "MSH6 IP, bandshift and MMR all WT-like", C_OLLILA)
put("msh2_msh6", "G322D", "binding", E3, "direct",
    "MMR-proficient with interaction intact; gnomAD ~1.9%", C_OLLILA)

# --------------------------------------------------------------- VHL-ElonginC
for ax in ("fold_complex", "binding"):
    put("vhl_elonginc", "L158Q", ax, E3, "direct",
        "Disrupts VCB complex assembly", C_STEBBINS)
    put("vhl_elonginc", "Y98H", ax, E3, "direct",
        "VCB assembly preserved while HIF binding is abolished "
        "(beta-domain HIF face, not the ElonginC face)", C_OHH)

# ------------------------------------------------------------------- KRAS-CRAF
# v8.0: the prior bases read "WT-like RAF1 binding" and cited Hunter 2015, but
# Hunter 2015 Table 2 ("Relative affinities for RAF kinase", nmol/L) reports
# WT 56 +/- 6, G12D 270 +/- 46, G12V 411 +/- 40, and the text states a "4.8-,
# 7.3-, and 6.2-fold decrease in affinity". The basis misread the source it
# named. Checked against the paper, not against a summary of it.
put("kras_craf", "G12D", "binding", E1, "direct",
    "RAF-RBD affinity decreased 4.8-fold (270 +/- 46 nmol/L against 56 +/- 6 "
    "for wild type). Measured magnitude 0.93 kcal/mol, below the 2.5 kcal/mol "
    "reference threshold", C_HUNTER)
put("kras_craf", "G12V", "binding", E1, "direct",
    "RAF-RBD affinity decreased 7.3-fold (411 +/- 40 nmol/L against 56 +/- 6 "
    "for wild type). Measured magnitude 1.18 kcal/mol, below the 2.5 kcal/mol "
    "reference threshold", C_HUNTER)
put("kras_craf", "Q61H", "binding", E3, "direct",
    "Preserved RAF binding; switch-II hydrolysis defect. Q61 is 11.1 A from RAF1",
    C_HUNTER)

# ------------------------------------------------------------------ hemoglobin
for ax, basis in [
    ("monomer", "Native beta-subunit fold; the lesion is intermolecular"),
    ("fold_complex", "Native tetramer; assembly unaffected"),
    ("binding", "HbS assembles normally -- pathogenicity arises from "
                "deoxy-state polymerization at a surface patch, not from "
                "alpha1-beta1 interface disruption"),
]:
    put("hemoglobin_dimer", "E6V", ax, E3, "direct", basis, C_EATON)

HB_DDG = {"W37Y": "+2.0", "W37A": "+5.0", "W37G": "+7.0", "W37E": "+9.0"}
for v, dd in HB_DDG.items():
    put("hemoglobin_tetramer", v, "binding", E1, "direct",
        f"Equilibrium dimer->tetramer assembly ddG = {dd} kcal/mol", C_KIGER)
    put("hemoglobin_tetramer", v, "fold_complex", E1, "coupled",
        "Same assembly-energy measurement; the fold-in-complex token is a "
        "coupled reading of the assembly ddG, not an independent "
        "subunit-stability measurement", C_KIGER)
# The source names this variant only as "hemoglobin Kansas" -- it never writes
# N102T -- and states the substitution as "threonine replaces asparagine in
# position 102 of the beta chain". It reports NO kcal value anywhere, only
# Kd(tetramer -> dimer): "about 1.0 to 1.4 x 10^-6 M" for normal oxyhemoglobin
# against "about 200 x 10^-6 M" for Hb Kansas. RT*ln of that ratio is 3.0
# kcal/mol per tetramer, or 1.5 per interface.
#
# The basis read "~1.5 kcal/mol", i.e. the per-interface reading, while the four
# beta-W37 rows beside it are per tetramer -- Kiger's own wording is "9
# kcal/tetramer" -- and all five sit in one calibration table under one axis and
# one system label. Harmonised to per tetramer, and the value is now DERIVED
# from these Kd figures in build_delta_calibration_points.py rather than typed.
put("hemoglobin_tetramer", "N102T", "binding", E1, "direct",
    "Hb Kansas: alpha1-beta2 H-bond loss. Source reports Kd(tetramer->dimer) "
    "~1.0-1.4 uM for normal oxyhemoglobin against ~200 uM, giving a reduced "
    "tetramer assembly energy of 3.0 kcal/mol per tetramer (1.5 per "
    "interface); stated per tetramer to match the beta-W37 series",
    C_BONAVENTURA)
put("hemoglobin_tetramer", "N102T", "fold_complex", E1, "coupled",
    "Coupled reading of the same assembly measurement", C_BONAVENTURA)

# ----------------------------------------------------------------- troponin I/C
put("troponin_ic", "R145G", "binding", E2, "direct",
    "EC50 1.4 uM vs 0.8 uM WT: reduced inhibitory potency, actin affinity "
    "unchanged -- interface retained", C_ELLIOTT)
put("troponin_ic", "R145Q", "binding", E2, "direct",
    "Same class as R145G: reduced inhibition, actin affinity unchanged",
    C_KIMURA)
for ax in AXES:
    put("troponin_ic", "P82S", ax, E4, "inferred",
        "AF ~3% African-American: benign by frequency", C_FRAZIER)

# ------------------------------------------------------------------ CaM-Cav1.2
for v in ("D96V", "N98S"):
    put("cam_cav12", v, "binding", E2, "off_axis",
        "Ca2+ affinity at EF-hand III is measured, but no IQ-motif binding "
        "measurement exists; the binding token is inferred from preserved "
        "target engagement", C_CROTTI)
put("cam_cav12", "F142L", "binding", E3, "coupled",
    "C-domain Ca2+ binding impaired while RyR2 inhibition is preserved -- "
    "engagement retained at a homologous IQ-type site", C_HWANG)

# ----------------------------------------------------------------- SMAD4-SMAD3
put("smad4_smad3", "D351H", "binding", E3, "direct",
    "Abolished phospho-Smad2/3 interaction by co-IP", C_DEBOSSCHER)
put("smad4_smad3", "R361C", "binding", E3, "direct",
    "Completely abolished R-Smad interaction", C_LANAUZE)
put("smad4_smad3", "I500T", "binding", E5, "inferred",
    "Post-translational mechanism (reduced ubiquitination -> protein "
    "stabilization); no interface measurement", C_LEGOFF)
# v8.0: I500V is no longer grouped with I500T. An interface measurement for
# this exact variant now exists -- Lindsay 2025 Fig 4D-E -- and it points the
# other way: the interaction is strengthened, not merely unmeasured. One of
# only two stabilizing commitments in the benchmark.
put("smad4_smad3", "I500V", "binding", E2, "direct",
    "GST pull-down shows increased recovery of SMAD3 with SMAD4 I500V, and "
    "TR-FRET shows a statistically significant increase in SMAD4-SMAD3 "
    "interaction. R361H and R361C negative controls show the opposite "
    "direction, so the assay resolves sign", C_LINDSAY)

# ------------------------------------------------------------------ CFH-C3b
for v, kd in [("R78G", "> 35 uM (binding effectively abolished)"),
              ("R53H", "~12 uM, comparable to WT"),
              ("I62V", "10-14 uM, WT-like")]:
    put("cfh_c3b", v, "monomer", E3, "direct",
        "Natively folded by CD / thermal denaturation", C_PECHTL)
    put("cfh_c3b", v, "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY,
        C_PECHTL)
    put("cfh_c3b", v, "binding", E1, "direct", f"SPR K_D {kd}", C_PECHTL)

# ---------------------------------------------------------------- VWF A1-GPIba
put("vwf_gpiba", "R1334Q", "monomer", E1, "direct",
    "Natively folded with modestly reduced Tm -- below the destabilization "
    "threshold", C_TISCHER)
put("vwf_gpiba", "R1334Q", "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY,
    C_TISCHER)
put("vwf_gpiba", "R1334Q", "binding", E1, "direct",
    "GPIba binding lost by SPR and shear-flow adhesion", C_TISCHER)
put("vwf_gpiba", "A1381T", "monomer", E1, "direct",
    "Natively folded, WT-like Tm", C_TISCHER)
put("vwf_gpiba", "A1381T", "fold_complex", E5, "inferred", NO_ASSEMBLY_ASSAY,
    C_TISCHER)
put("vwf_gpiba", "A1381T", "binding", E1, "direct",
    "Binding-competent by SPR", C_TISCHER)


def brct_assignments() -> dict:
    """BRCT monomer axes: measured GdmCl unfolding energies, read from file.

    Never transcribed -- the measured value is pulled from the concordance
    file so a re-curation propagates.
    """
    out = {}
    b = pd.read_csv(BRCT)
    for _, r in b.iterrows():
        val = r["measured_ddG_UF_kcal_mol"]
        out[("brca1_brct", r["variant"], "monomer")] = (
            E1, "direct",
            f"GdmCl equilibrium unfolding: measured ddG_U-F = {val:+.2f} kcal/mol",
            C_ROWLING,
        )
    return out


def main() -> int:
    df = pd.read_csv(CANON)
    ASSIGN.update(brct_assignments())

    rows = []
    for _, r in df.iterrows():
        for axis, col in AXES.items():
            token = r[col]
            if token not in COMMITTED:
                continue
            key = (r["system"], r["variant"], axis)
            if key not in ASSIGN:
                print(f"UNCURATED AXIS: {key}", file=sys.stderr)
                return 1
            etype, direct, basis, cite = ASSIGN[key]
            rows.append(dict(
                system=r["system"], variant=r["variant"], axis=axis,
                expected_token=token, evidence_type=etype,
                evidence_directness=direct, evidence_basis=basis,
                evidence_citation=cite,
            ))
    led = pd.DataFrame(rows).sort_values(["system", "variant", "axis"])

    extra = set(ASSIGN) - {(r.system, r.variant, r.axis)
                           for r in led.itertuples()}
    # An ASSIGN entry for an axis the canonical marks `unknown` is expected when that
    # axis was deliberately withdrawn (v7.6/v7.7): the curated evidence record is kept
    # on purpose as the audit trail for why the axis is ungraded. Only unexplained
    # entries are a curation error.
    withdrawn = set()
    wpath = REPO / "reference_outputs" / "COMAVI_evidence_ledger_withdrawn_v76.csv"
    if wpath.exists():
        w = pd.read_csv(wpath)
        withdrawn = {(r.system, r.variant, r.axis) for r in w.itertuples()}
    unexplained = sorted(extra - withdrawn)
    if extra & withdrawn:
        print(f"note: {len(extra & withdrawn)} ASSIGN entries retained for withdrawn "
              f"axes (audit trail); not committed", file=sys.stderr)
    if unexplained:
        print(f"ASSIGN entries for non-committed axes: {unexplained}",
              file=sys.stderr)
        return 1

    led.to_csv(OUT_CSV, index=False)

    by_type = led.evidence_type.value_counts().to_dict()
    by_dir = led.evidence_directness.value_counts().to_dict()
    xt = (led.groupby(["evidence_type", "evidence_directness"]).size()
          .unstack(fill_value=0))
    quant = led.evidence_type.isin([E1, E2])
    axis_ok = led.evidence_directness == "direct"
    summary = dict(
        n_committed_axes=int(len(led)),
        n_variants=int(led.variant.nunique()),
        n_systems=int(led.system.nunique()),
        by_evidence_type=by_type,
        by_directness=by_dir,
        by_axis={a: int(n) for a, n in led.axis.value_counts().items()},
        n_energetic_direct=int(((led.evidence_type == E1) & axis_ok).sum()),
        frac_energetic=round(float((led.evidence_type == E1).mean()), 4),
        frac_quantitative=round(float(quant.mean()), 4),
        frac_direct=round(float(axis_ok.mean()), 4),
        frac_inferred=round(float((led.evidence_type == E5).mean()), 4),
        type_by_axis={
            a: g.evidence_type.value_counts().to_dict()
            for a, g in led.groupby("axis")
        },
        note=("Denominator is COMMITTED axes (expected_ddg_* in "
              "{destab,neutral,stab}); unknown/not_applicable axes are not "
              "scored and carry no ground truth to characterize."),
    )
    OUT_JSON.write_text(json.dumps(summary, indent=1))

    print(f"wrote {OUT_CSV.relative_to(REPO)}  ({len(led)} committed axes)")
    print(f"wrote {OUT_JSON.relative_to(REPO)}")
    print("\ntype x directness:")
    print(xt.to_string())
    print("\nby axis:")
    print(led.groupby(["axis", "evidence_type"]).size().unstack(fill_value=0)
          .to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
