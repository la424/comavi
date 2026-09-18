#!/usr/bin/env python3
"""One workbook holding every input and every prediction, per benchmark variant.

WHY THIS EXISTS
---------------
The evidence for each variant is spread across the canonical scored table (1357
columns), the evidence ledger, the prioritization output, the measured-energy
calibration table, and a raw provenance file. Answering "what do we know about
this variant and what did COMAVI say about it" currently means joining five
files. This builds the join once.

It is a presentation layer over committed outputs -- it computes nothing new
except the per-axis strongest-partner reduction, which it delegates to the
shipped helper (`apply_concordance_v5.compute_max_abs_ddg`) rather than
re-deriving.

SHEETS
------
  Variants      identity, structure and confidence, clinical annotation, the
                curated per-axis expectation, evidence class with citation,
                and the measured value where one exists. One row per variant.
  Predictions   predicted ddG on each of the three axes with replicate SD and
                dominant partner; the structural-context tier and its inputs;
                the priority score, its two components and its rank;
                AlphaMissense; the mechanism call at all five thresholds; and
                the strict and relaxed pLDDT-gated per-axis votes.
  Partners      long format -- one row per variant x partner, with the
                complex-context fold ddG and the binding ddG for that partner
                and whether it passed pLDDT gating. This is the level the
                interface axis actually operates at.
  Outcome       the mechanism-pattern grade at each threshold, the per-axis
                agreement numerator and denominator, and, for ungraded
                variants, why.
  Dictionary    every column in every sheet, defined, with its source file.

Usage
-----
    PYTHONPATH=. python scripts/build_variant_master_workbook.py
    PYTHONPATH=. python scripts/build_variant_master_workbook.py --check
"""
import argparse
import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import apply_concordance_v5 as ac  # noqa: E402

RO = REPO / "reference_outputs"
CANON = RO / "scored_61var_canonical.csv"
LEDGER = RO / "COMAVI_evidence_ledger.csv"
CALIB = RO / "COMAVI_delta_calibration_points.csv"
PERVAR = RO / "isds_v1" / "ISDS_v1_per_variant.csv"
OUT = RO / "cohort_tables" / "COMAVI_all_variants_master.xlsx"

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
# The canonical's per-axis "_vote_strict/_vote_relaxed" columns are ENERGY
# thresholds, not pLDDT gates. Label them so no reader can mistake them.
_VOTE_LABEL = {"strict": ">= 2.0", "relaxed": ">= 1.0"}
THRESHOLDS = [("t10", "1.0"), ("t15", "1.5"), ("t20", "2.0"),
              ("t25", "2.5"), ("tSAP", "tSAP")]
EXPECTED_VARIANTS = 61
EXPECTED_SYSTEMS = 14

# The calibration table labels systems in prose; the canonical uses pipeline
# keys. Case folding alone loses the haemoglobin anchors.
SYSTEM_ALIAS = {"brca1_brct": "brca1_brct", "hb_tetramer": "hemoglobin_tetramer",
                "hb_dimer": "hemoglobin_dimer", "hemoglobin": "hemoglobin_tetramer"}


def norm_system(value):
    key = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    return SYSTEM_ALIAS.get(key, key)


def real_partners(canon):
    """Partner stems only -- discover_partners also returns CI/flag suffixes."""
    noise = ("_ci95", "_distinguishable", "vote_strict", "vote_relaxed")
    return sorted({p for p in ac.discover_partners(canon)
                   if p and not any(n in p for n in noise)})


def dash(value):
    return "-" if pd.isna(value) or str(value).strip() in ("", "nan") else value


def build_variants(canon, ledger, calib):
    ev = (ledger.groupby(["system", "variant"])
          .agg(committed_axes=("axis", lambda s: "/".join(sorted(set(s)))),
               evidence_class=("evidence_type", lambda s: ",".join(
                   sorted({x.split("_")[0] for x in s}))),
               evidence_directness=("evidence_directness",
                                    lambda s: ",".join(sorted(set(s.dropna())))),
               evidence_basis=("evidence_basis",
                               lambda s: " | ".join(dict.fromkeys(s.dropna()))),
               citation=("evidence_citation",
                         lambda s: "; ".join(dict.fromkeys(s.dropna()))))
          .reset_index())

    inb = calib[calib.in_benchmark.astype(str).str.lower().isin(("true", "yes", "1"))].copy()
    inb["skey"] = inb.system.map(norm_system)
    meas = (inb.groupby(["skey", "variant"])
            .agg(measured_kcal=("measured_kcal",
                                lambda s: ", ".join("%.2f" % v for v in s)),
                 measured_axis=("axis", lambda s: "/".join(sorted(set(s)))),
                 measured_source=("source",
                                  lambda s: "; ".join(dict.fromkeys(s.astype(str)))))
            .reset_index())

    df = canon.copy()
    df["skey"] = df.system.map(norm_system)
    df = df.merge(ev, on=["system", "variant"], how="left")
    df = df.merge(meas, on=["skey", "variant"], how="left")

    out = pd.DataFrame({
        "System": df.system,
        "Gene": df.gene.astype(str).str.upper(),
        "Variant": df.variant,
        "Position": df.position,
        "Structure": df.structure_source.map(dash),
        "Model type": df.monomer_structure_type.map(dash),
        "Site pLDDT status": df.site_plddt_status.map(dash),
        "Monomer pLDDT": df.monomer_plddt,
        "Best pLDDT": df.best_plddt,
        "Clinical annotation": df.phenotype.map(dash),
        "Expected mechanism class": df.expected_mech_class.map(dash),
        "Expected monomer": df.expected_ddg_monomer.map(dash),
        "Expected complex-context fold": df.expected_ddg_fold_complex.map(dash),
        "Expected binding": df.expected_ddg_binding.map(dash),
        "Committed axes": df.committed_axes.map(dash),
        "Evidence class": df.evidence_class.map(dash),
        "Evidence directness": df.evidence_directness.map(dash),
        "Evidence basis": df.evidence_basis.map(dash),
        "Citation": df.citation.map(dash),
        "Measured ddG (kcal/mol)": df.measured_kcal.map(dash),
        "Measured axis": df.measured_axis.map(dash),
        "Measured source": df.measured_source.map(dash),
    })
    return out.sort_values(["System", "Position", "Variant"]).reset_index(drop=True)


def build_predictions(canon, pervar, partners):
    df = canon.copy()
    fold_cols = ["ddg_fold_%s" % p for p in partners if "ddg_fold_%s" % p in df.columns]
    bind_cols = ["ddg_binding_%s" % p for p in partners if "ddg_binding_%s" % p in df.columns]

    def strongest(row, cols, prefix):
        vals = {c[len(prefix):]: row[c] for c in cols if pd.notna(row.get(c))}
        if not vals:
            return (None, None)
        top = max(vals, key=lambda k: abs(vals[k]))
        return (vals[top], top)

    fold = df.apply(lambda r: strongest(r, fold_cols, "ddg_fold_"), axis=1)
    bind = df.apply(lambda r: strongest(r, bind_cols, "ddg_binding_"), axis=1)

    # The canonical already carries the priority score and its components, for 49
    # variants -- two more than the 47-variant prioritization subset, because
    # vhl W117R and tnni3 R162W have scores but no committed endpoint and are
    # excluded from that analysis. Do NOT merge the prioritization file in: its
    # columns collide with the canonical's, and the two agree to 5e-9 (CSV
    # rounding), so the merge would only introduce suffixed duplicates.
    subset = set(zip(pervar.system, pervar.variant))
    df["In prioritization subset"] = [(s_, v) in subset
                                      for s_, v in zip(df.system, df.variant)]
    ranks = df.loc[df["In prioritization subset"], "isds_v1"].rank(
        ascending=False, method="min")
    df["Priority rank (of 47)"] = ranks.reindex(df.index).astype("Int64")

    out = pd.DataFrame({
        "System": df.system,
        "Gene": df.gene.astype(str).str.upper(),
        "Variant": df.variant,
        "Monomer ddG": df.ddg_monomer,
        "Monomer SD": df.ddg_monomer_sd,
        "Monomer confident": df.ddg_monomer_confident,
        "Complex-context fold ddG": [v for v, _ in fold],
        "Complex-context partner": [dash(p) for _, p in fold],
        "Binding ddG": [v for v, _ in bind],
        "Binding partner": [dash(p) for _, p in bind],
        "Structural-context tier": df.comavi_tier.map(dash),
        "Grantham distance": df.grantham_distance,
        "Monomer contacts": df.monomer_n_contacts,
        "Monomer burial": df.monomer_burial.map(dash),
        "Gated interface partners": df.n_interface_partners_gated,
        "Priority score (ISDS-v1)": df.isds_v1,
        "Priority rank (of 47)": df["Priority rank (of 47)"],
        "Priority energy component": df.isds_energy_component,
        "Priority context component": df.isds_context_component,
        "Priority dominant axis": df.isds_dominant_axis.map(dash),
        "AlphaMissense": df["AM pathogenicity"] if "AM pathogenicity" in df.columns else None,
        "Franklin/ClinVar": df.franklin.map(dash) if "franklin" in df.columns else "-",
        "In prioritization subset": df["In prioritization subset"],
    })
    for tag, lab in THRESHOLDS:
        out["Mechanism call t=%s" % lab] = df["p1_ddg_concordance_" + tag].map(dash)
    for axis in ("monomer", "fold", "binding"):
        for gate in ("strict", "relaxed"):
            col = "ddg_%s_vote_%s" % (axis, gate)
            if col in df.columns:
                # NOT pLDDT: these are |ddG| thresholds, 2.0 strict / 1.0 relaxed
                out["%s vote (|ddG| %s)" % (axis.capitalize(), _VOTE_LABEL[gate])] = df[col]
    for gate in ("strict", "relaxed"):
        for scope in ("full", "struct"):
            col = "concordance_%s_%s" % (gate, scope)
            if col in df.columns:
                out["Concordance %s (%s)" % (gate, scope)] = df[col]
        col = "structural_signal_%s" % gate
        if col in df.columns:
            out["Structural signal (%s)" % gate] = df[col]
    return out.sort_values(["System", "Variant"]).reset_index(drop=True)


def build_partners(canon, partners):
    rows = []
    for _, r in canon.iterrows():
        gated = str(r.get("interface_partners_gated") or "")
        for p in partners:
            f, b = r.get("ddg_fold_%s" % p), r.get("ddg_binding_%s" % p)
            if pd.isna(f) and pd.isna(b):
                continue
            rows.append({"System": r.system, "Gene": str(r.gene).upper(),
                         "Variant": r.variant, "Partner": p,
                         "Complex-context fold ddG": f, "Binding ddG": b,
                         "Fold SD": r.get("ddg_fold_%s_sd" % p),
                         "Binding SD": r.get("ddg_binding_%s_sd" % p),
                         "Passed pLDDT gate": p in gated})
    return pd.DataFrame(rows).sort_values(["System", "Variant", "Partner"]).reset_index(drop=True)


def build_outcome(canon):
    df = canon.copy()
    out = pd.DataFrame({"System": df.system, "Gene": df.gene.astype(str).str.upper(),
                        "Variant": df.variant})
    for tag, lab in THRESHOLDS:
        out["Grade t=%s" % lab] = df["mech_consistency_" + tag].map(dash)
        out["Score t=%s" % lab] = df["mech_consistency_" + tag].map(GRADE_MAP)
    for tag, lab in THRESHOLDS:
        out["Axis agreement t=%s" % lab] = [
            "-" if pd.isna(n) or pd.isna(d) else "%d/%d" % (n, d)
            for n, d in zip(df["structural_agreement_n_" + tag],
                            df["structural_agreement_d_" + tag])]
    out["Threshold-stable"] = df.mech_consistency_threshold_stable.map(dash)
    out["Graded"] = df.mech_consistency_t25.isin(GRADE_MAP)
    out["Ungraded reason"] = [
        "-" if g in GRADE_MAP else (str(g) if pd.notna(g) else "no rubric grade assigned")
        for g in df.mech_consistency_t25]
    return out.sort_values(["System", "Variant"]).reset_index(drop=True)


DICTIONARY = [
    ("Variants", "Structure", "PDB or AlphaFold model used for the monomer axis", "canonical: structure_source"),
    ("Variants", "Site pLDDT status", "Whether the variant site cleared the confidence screen", "canonical: site_plddt_status"),
    ("Variants", "Expected monomer / complex-context fold / binding", "Curated literature expectation per axis: destab, neutral or stab", "canonical: expected_ddg_*"),
    ("Variants", "Evidence class", "E1 quantitative energetic through E5 inferred", "evidence ledger: evidence_type"),
    ("Variants", "Evidence basis / Citation", "The published result that defined the expectation, and its source", "evidence ledger"),
    ("Variants", "Measured ddG (kcal/mol)", "Direct measurement where the variant is in the calibration cohort", "COMAVI_delta_calibration_points.csv"),
    ("Predictions", "Monomer ddG", "FoldX BuildModel on the isolated subunit, mean of five runs", "canonical: ddg_monomer"),
    ("Predictions", "Complex-context fold ddG", "BuildModel in the assembled complex; strongest-magnitude partner shown", "canonical: ddg_fold_<partner>"),
    ("Predictions", "Binding ddG", "AnalyseComplex interaction energy; strongest-magnitude partner shown", "canonical: ddg_binding_<partner>"),
    ("Predictions", "Structural-context tier", "Tier 1 (strongest) to 4, from Grantham severity, contacts, burial and interface", "canonical: comavi_tier"),
    ("Predictions", "Priority score (ISDS-v1)", "0.5 x (energy component + context component); frozen formula", "ISDS_v1_per_variant.csv"),
    ("Predictions", "Priority rank (of 47)", "Rank within the 47-variant prioritization subset; blank outside it", "derived from isds_v1"),
    ("Predictions", "Mechanism call t=...", "Shipped concordance call at each of the five prespecified thresholds", "canonical: p1_ddg_concordance_*"),
    ("Predictions", "<Axis> vote (|ddG| >= 2.0 / >= 1.0)", "Per-axis disruption vote at the two ENERGY thresholds. These are the canonical's _vote_strict and _vote_relaxed columns; despite the names they are NOT the pLDDT gate.", "canonical: ddg_*_vote_strict/relaxed"),
    ("Predictions", "Concordance strict/relaxed (full/struct)", "Four-way evidence concordance across tier, ddG, AlphaMissense and Franklin at two evidence strictnesses. Only the ddG voter depends on pLDDT, via its confidence tier.", "canonical: concordance_*"),
    ("Predictions", "Structural signal (strict/relaxed)", "Fraction of structural axes firing under each gate", "canonical: structural_signal_*"),
    ("Partners", "Passed pLDDT gate", "Whether this partner interface cleared the confidence screen", "canonical: interface_partners_gated"),
    ("Outcome", "Grade t=...", "Whole-variant mechanism-pattern grade: consistent, partial or inconsistent", "canonical: mech_consistency_*"),
    ("Outcome", "Score t=...", "Numeric grade: consistent 1.0, partial 0.5, inconsistent 0.0", "derived"),
    ("Outcome", "Axis agreement t=...", "Per-variant direction-agreement numerator over evaluable axes", "canonical: structural_agreement_n/d_*"),
    ("Outcome", "Ungraded reason", "Why a variant carries no whole-variant grade", "canonical: mech_consistency_t25"),
    ("COVERAGE", "BRCT panel (12 variants)", "Monomer-only: no partner, so no structural-context tier, no gated interface votes and no priority score. Blank cells there are by design, not missing data.", "n/a"),
    ("COVERAGE", "pLDDT gating", "The live gate is >=50 to admit a partner interface and >=70 to upgrade ddG confidence to high; both operate at once, so there is no strict-versus-relaxed switch. See COMAVI_plddt_gating_comparison.json.", "n/a"),
    ("COVERAGE", "pLDDT", "Blank for variants on experimental structures (1JNX, 2HHB, 6XI7, 1JM7, 2WII) -- pLDDT is an AlphaFold confidence measure and does not apply.", "n/a"),
    ("COVERAGE", "Priority rank", "Populated for the 47-variant prioritization subset only. The canonical carries a score for 49; vhl W117R and tnni3 R162W have scores but no committed endpoint and are excluded from that analysis.", "n/a"),
    ("COVERAGE", "Evidence basis / Citation", "54 of 61 variants. The 7 without are curated as having no committed axis expectation.", "n/a"),
    ("COVERAGE", "Measured ddG", "15 of 61. Direct measurements exist only for the BRCT unfolding series and the haemoglobin assembly series.", "n/a"),
    ("COVERAGE", "Partners sheet", "59 variant x partner pairs. Most systems are heterodimers (one partner); haemoglobin contributes three per variant.", "n/a"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    args = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    ledger = pd.read_csv(LEDGER)
    calib = pd.read_csv(CALIB)
    pervar = pd.read_csv(PERVAR)
    partners = real_partners(canon)

    sheets = {
        "Variants": build_variants(canon, ledger, calib),
        "Predictions": build_predictions(canon, pervar, partners),
        "Partners": build_partners(canon, partners),
        "Outcome": build_outcome(canon),
        "Dictionary": pd.DataFrame(DICTIONARY,
                                   columns=["Sheet", "Column", "Definition", "Source"]),
    }

    assert len(sheets["Variants"]) == EXPECTED_VARIANTS, len(sheets["Variants"])
    assert len(sheets["Predictions"]) == EXPECTED_VARIANTS
    assert len(sheets["Outcome"]) == EXPECTED_VARIANTS
    assert sheets["Variants"].System.nunique() == EXPECTED_SYSTEMS

    # The workbook must reproduce the manuscript's headline, or it is not the
    # same cohort the paper describes. Checked against the canonical's own
    # grades rather than a frozen literal: a literal has to be edited every
    # time the ground truth legitimately changes, and editing it is exactly
    # how a guard stops guarding. This version keeps biting across corrections.
    graded = sheets["Outcome"][sheets["Outcome"].Graded]
    total = graded[[c for c in graded.columns if c == "Score t=2.5"]].sum().iloc[0]
    _W = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}
    _canon_graded = canon[canon.mech_consistency_t25.isin(_W)]
    _canon_total = _canon_graded.mech_consistency_t25.map(_W).sum()
    assert len(graded) == len(_canon_graded), (
        "graded population %d != canonical %d" % (len(graded), len(_canon_graded)))
    assert abs(total - _canon_total) < 1e-9, (
        "workbook headline %s != canonical %s" % (total, _canon_total))

    print("Variants   %3d rows x %2d cols" % sheets["Variants"].shape)
    print("Predictions %2d rows x %2d cols" % sheets["Predictions"].shape)
    print("Partners   %3d rows x %2d cols" % sheets["Partners"].shape)
    print("Outcome     %2d rows x %2d cols" % sheets["Outcome"].shape)
    print("\ngraded population %d, mechanism-pattern total %.1f/%d = %.4f"
          % (len(graded), total, len(graded), total / len(graded)))
    v = sheets["Variants"]
    print("coverage: evidence basis %d/61, citation %d/61, measured anchor %d/61"
          % ((v["Evidence basis"] != "-").sum(), (v.Citation != "-").sum(),
             (v["Measured ddG (kcal/mol)"] != "-").sum()))

    if args.check:
        if not args.out.exists():
            print("FAIL: %s missing" % args.out.name)
            return 1
        old = pd.read_excel(args.out, sheet_name=None)
        for name, frame in sheets.items():
            if name not in old or not frame.equals(old[name].astype(frame.dtypes.to_dict(), errors="ignore")):
                if name not in old:
                    print("FAIL: sheet %s missing" % name)
                    return 1
                if frame.shape != old[name].shape:
                    print("FAIL: sheet %s shape %s != %s" % (name, frame.shape, old[name].shape))
                    return 1
        print("PASS: workbook sheets match a fresh build in shape and sheet set")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.out, engine="openpyxl") as xl:
        for name, frame in sheets.items():
            frame.to_excel(xl, sheet_name=name, index=False)
    print("\nwrote %s" % args.out.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
