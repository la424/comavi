#!/usr/bin/env python3
"""Build the benchmark cohort tables the reviewers asked for.

WHY THIS EXISTS
---------------
The manuscript shipped fifteen tables and none of them listed a variant. Two
lab readers independently could not tell what the benchmark was: "one mutation
in one gene, multiple mutations in one gene, or different mutations in
different genes?", "a set of mutations that was previously specified and used
for all the benchmarking tests?", "why this mutation? anything important about
it?". Six of one reviewer's fifteen comments reduce to that single absence.

Two tables answer it, and both are derivable from committed data -- nothing
here is authored:

  Table 1  (main text)   one row per protein system: what was modeled, from
                         what structure, how many variants, how many gradeable.
  Table S1 (supplement)  one row per variant: clinical annotation, the expected
                         mechanism (ground truth), the structural-context tier,
                         COMAVI's graded agreement, and -- the part that answers
                         "why this variant" -- the one-line evidence basis and
                         its literature citation, both already curated in
                         reference_outputs/COMAVI_evidence_ledger.csv.

SOURCES
    reference_outputs/scored_61var_canonical.csv     cohort, labels, tiers
    reference_outputs/COMAVI_evidence_ledger.csv     basis + citation per axis
    reference_outputs/COMAVI_delta_calibration_points.csv   measured anchors

USAGE
    python scripts/build_variant_tables.py                 # write both CSVs
    python scripts/build_variant_tables.py --check          # exit 1 on drift
    python scripts/build_variant_tables.py --print-only
"""

import argparse
import io
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import apply_concordance_v5 as _AC

REPO = Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
LEDGER = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"
CALIB = REPO / "reference_outputs" / "COMAVI_delta_calibration_points.csv"
BRCT = REPO / "supplement" / "brct" / "brct_foldx_concordance.csv"
OUT_DIR = REPO / "reference_outputs" / "cohort_tables"

EXPECTED_SYSTEMS = 14
EXPECTED_VARIANTS = 61

# The committed threshold set, in order. Sourced from the pipeline's own spec so
# a threshold added there appears here without editing this file.
THRESHOLD_TAGS = [tag for tag, _ in _AC.THRESHOLD_SPECS]

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}

# Display names. The canonical's snake_case keys are pipeline identifiers, not
# something a reader should have to decode.
SYSTEM_LABEL = {
    "brca1_bard1": "BRCA1-BARD1 RING heterodimer",
    "brca1_brct": "BRCA1 BRCT tandem domain (fold only)",
    "cam_cav12": "Calmodulin-Ca(v)1.2 IQ domain",
    "cfh_c3b": "Complement factor H-C3b",
    "hemoglobin_dimer": "Hemoglobin alpha-beta dimer",
    "hemoglobin_tetramer": "Hemoglobin tetramer",
    "kras_craf": "KRAS-CRAF RBD",
    "mlh1_pms2": "MLH1-PMS2 C-terminal heterodimer",
    "msh2_msh6": "MSH2-MSH6 heterodimer",
    "pi3k": "PI3K p110alpha-p85alpha",
    "smad4_smad3": "SMAD4-SMAD3 MH2 complex",
    "troponin_ic": "Troponin I-troponin C",
    "vhl_elonginc": "VHL-elongin C",
    "vwf_gpiba": "VWF A1-GPIb alpha",
}

MECH_LABEL = {
    "fold_mechanism": "Fold destabilization",
    "ppi_destab_mechanism": "Interface disruption",
    "mixed_structural": "Fold and interface",
    "structurally_silent": "No modeled lesion",
    "structurally_uncommitted": "Not committed",
    "interface_uncommitted_magnitude": "Interface, magnitude uncommitted",
}

PHENO_LABEL = {
    "pathogenic": "Pathogenic",
    "pathogenic_gof": "Pathogenic (gain of function)",
    "benign": "Benign",
    "vus": "Uncertain significance",
}


# The calibration table names systems in prose ('BRCA1 BRCT', 'Hb tetramer')
# while the canonical uses pipeline keys. Case folding alone recovers only
# 10 of the 15 in-benchmark anchors; these aliases recover the rest. Without
# them the haemoglobin binding anchors drop silently -- a join that loses rows
# and reports success.
SYSTEM_ALIAS = {
    "hb_tetramer": "hemoglobin_tetramer",
    "hb_dimer": "hemoglobin_dimer",
    "hemoglobin": "hemoglobin_tetramer",
}


def norm_system(value):
    """Canonical system key for any of the labels used across the inputs."""
    key = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    return SYSTEM_ALIAS.get(key, key)


def pretty(value, mapping):
    if pd.isna(value):
        return "-"
    key = str(value)
    return mapping.get(key, key.replace("_", " "))


def load():
    canonical = pd.read_csv(CANON, low_memory=False)
    ledger = pd.read_csv(LEDGER)
    calib = pd.read_csv(CALIB)
    return canonical, ledger, calib


def evidence_per_variant(ledger):
    """One row per variant: axes committed, evidence classes, basis, citation."""
    def uniq(series, sep):
        return sep.join(dict.fromkeys(str(x) for x in series.dropna()))

    grouped = ledger.groupby(["system", "variant"]).agg(
        committed_axes=("axis", lambda s: uniq(s.sort_values(), ", ")),
        evidence_class=("evidence_type",
                        lambda s: uniq(s.map(lambda x: str(x).split("_")[0]).sort_values(), "/")),
        evidence_directness=("evidence_directness", lambda s: uniq(s.sort_values(), "/")),
        evidence_basis=("evidence_basis", lambda s: uniq(s, "; ")),
        citation=("evidence_citation", lambda s: uniq(s, "; ")),
    ).reset_index()
    return grouped


def measured_per_variant(calib):
    """Measured energetic anchors for variants that have one."""
    frame = calib.copy()
    flag = frame["in_benchmark"].astype(str).str.lower().isin(("true", "yes", "1"))
    frame = frame[flag]
    frame["system"] = frame["system"].map(norm_system)
    return frame.groupby(["system", "variant"]).agg(
        measured_kcal=("measured_kcal", lambda s: ", ".join(f"{v:.2f}" for v in s)),
        measured_axis=("axis", lambda s: ", ".join(dict.fromkeys(str(x) for x in s))),
    ).reset_index()


def brct_measured(frame):
    """Measured GdmCl unfolding energies for the BRCT panel.

    The kcal/mol calibration table carries only the BRCT rows that enter the
    predicted-vs-measured comparison, which excludes the fold-intact /
    function-lost class by design -- that class is not a stability test. Those
    variants nevertheless HAVE a measured unfolding energy, and reporting it is
    the whole point of a per-variant table, so it is read from the panel's own
    source here rather than left blank.
    """
    if not BRCT.is_file():
        return pd.DataFrame(columns=["system_key", "variant", "brct_measured",
                                     "brct_role"])
    src = pd.read_csv(BRCT)
    have = set(zip(frame["system_key"], frame["variant"]))
    out = []
    for _, r in src.iterrows():
        if ("brca1_brct", r["variant"]) not in have:
            continue
        v = r.get("measured_ddG_UF_kcal_mol")
        if pd.isna(v):
            continue
        out.append({"system_key": "brca1_brct", "variant": r["variant"],
                    "brct_measured": "%.2f" % float(v),
                    "brct_role": str(r.get("role_in_cohort", "-"))})
    return pd.DataFrame(out)


def predicted_per_axis(canonical):
    """Strongest predicted energy on each axis, per variant.

    Delegates to the shipped helpers rather than re-deriving: axis_columns and
    strongest are the same functions that decide which axis fires in the
    detection and attribution metrics, so a value shown here cannot disagree
    with the call the paper reports for that variant.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import build_detection_vs_attribution as dva

    fold_cols = dva.axis_columns(canonical, "fold")
    bind_cols = dva.axis_columns(canonical, "binding")
    rows = []
    for _, r in canonical.iterrows():
        mono = r.get("ddg_monomer")
        fold = dva.strongest(r, fold_cols, "ddg_fold_")
        bind = dva.strongest(r, bind_cols, "ddg_binding_")
        rows.append({
            "variant_key": (norm_system(r["system"]), r["variant"]),
            "pred_monomer": "-" if pd.isna(mono) else "%.2f" % float(mono),
            "pred_fold": "-" if fold is None or pd.isna(fold) else "%.2f" % float(fold),
            "pred_binding": "-" if bind is None or pd.isna(bind) else "%.2f" % float(bind),
        })
    return pd.DataFrame(rows)


def agreement_across_thresholds(canonical):
    """The graded call at every committed threshold, not only the reference.

    A single reference-threshold grade hides whether a variant is graded the
    same way everywhere or flips; the sweep in Section 3.5 is an aggregate over
    exactly these per-variant columns.
    """
    cols = ["mech_consistency_%s" % tag for tag in THRESHOLD_TAGS]
    missing = [c for c in cols if c not in canonical.columns]
    assert not missing, "canonical lacks graded columns: %s" % missing
    short = {"consistent": "C", "partial": "P", "inconsistent": "X"}
    out = []
    for _, r in canonical.iterrows():
        marks = [short.get(str(r[c]), "-") for c in cols]
        graded = any(m != "-" for m in marks)
        out.append({
            "variant_key": (norm_system(r["system"]), r["variant"]),
            "agreement_sweep": "/".join(marks),
            # An ungraded variant has no call to be stable, so it gets a dash
            # rather than "yes" -- reporting agreement where there is no
            # judgement would inflate the stable count.
            "agreement_stable": ("yes" if len(set(marks)) == 1 else "no")
                                if graded else "-",
        })
    return pd.DataFrame(out)


def whole_variant_expectation(canonical):
    """Basis and citation for variants whose expectation is not per-axis.

    Five graded variants commit no axis at all: their curated expectation is a
    WHOLE-VARIANT one -- pathogenic by a functional mechanism, therefore no
    structural lesion should be detected -- carried by evidence_axes and the
    curation note. The per-axis evidence ledger is keyed on committed axis
    tokens, so it holds no row for them, and joining on it alone leaves their
    basis and citation blank. That blank reads as missing curation when the
    curation exists; it is just recorded at a different granularity.
    """
    rows = []
    for _, r in canonical.iterrows():
        note = str(r.get("notes") or "").strip()
        if note.lower() in ("", "nan"):
            continue
        rows.append({"variant_key": (norm_system(r["system"]), r["variant"]),
                     "note_basis": re.sub(r"\s+", " ", note),
                     "evidence_axes": str(r.get("evidence_axes") or "-")})
    return pd.DataFrame(rows)


def variant_table(canonical, ledger, calib):
    frame = canonical.copy()
    frame["system_key"] = frame["system"].map(norm_system)

    evidence = evidence_per_variant(ledger)
    evidence["system_key"] = evidence["system"].map(norm_system)
    measured = measured_per_variant(calib)
    measured = measured.rename(columns={"system": "system_key"})

    merged = frame.merge(evidence.drop(columns=["system"]),
                         on=["system_key", "variant"], how="left")
    merged = merged.merge(measured, on=["system_key", "variant"], how="left")

    brct = brct_measured(frame)
    if len(brct):
        merged = merged.merge(brct, on=["system_key", "variant"], how="left")
    else:
        merged["brct_measured"] = pd.NA
        merged["brct_role"] = pd.NA

    key = list(zip(merged["system_key"], merged["variant"]))
    pred = predicted_per_axis(canonical).set_index("variant_key")
    agree = agreement_across_thresholds(canonical).set_index("variant_key")
    notes = whole_variant_expectation(canonical).set_index("variant_key")

    def lookup(table, col, default="-"):
        return [table[col].get(k, default) for k in key]

    # A measured energy the kcal/mol calibration table does not carry, because
    # that table holds only the rows entering the predicted-vs-measured
    # comparison. Fall back to the panel's own source before writing a dash.
    meas = []
    for i, k in enumerate(key):
        v = merged["measured_kcal"].iloc[i]
        if isinstance(v, str) and v.strip() and v != "-":
            meas.append(v)
            continue
        b = merged["brct_measured"].iloc[i]
        meas.append(b if isinstance(b, str) and b.strip() else "-")

    # Basis and citation for the whole-variant expectations the per-axis ledger
    # cannot hold. Marked so a reader can see the granularity differs.
    basis = []
    cites = []
    for i, k in enumerate(key):
        b = merged["evidence_basis"].iloc[i]
        c = merged["citation"].iloc[i]
        if isinstance(b, str) and b.strip() and b != "-":
            basis.append(b)
            cites.append(c if isinstance(c, str) and c.strip() else "-")
            continue
        nb = notes["note_basis"].get(k)
        if nb:
            basis.append("Whole-variant expectation (no axis committed): %s" % nb)
            cites.append("see curation note")
        else:
            basis.append("-")
            cites.append(c if isinstance(c, str) and c.strip() else "-")

    out = pd.DataFrame({
        "System": merged["system_key"].map(lambda k: SYSTEM_LABEL.get(k, k)),
        "Gene": merged["gene"].fillna("-"),
        "Variant": merged["variant"],
        "Clinical annotation": merged["phenotype"].map(lambda v: pretty(v, PHENO_LABEL)),
        "Expected mechanism": merged["expected_mech_class"].map(lambda v: pretty(v, MECH_LABEL)),
        "Structural tier": merged["comavi_tier"].fillna("-"),
        "COMAVI agreement (t=2.5)": merged["mech_consistency_t25"].fillna("ungraded"),
        "Agreement across thresholds": lookup(agree, "agreement_sweep"),
        "Same call at every threshold": lookup(agree, "agreement_stable"),
        "Grantham": merged["grantham_distance"],
        "Gated interface partners": merged["n_interface_partners_gated"],
        "Committed axes": merged["committed_axes"].fillna("none committed"),
        "Evidence class": merged["evidence_class"].fillna("-"),
        "Evidence directness": merged["evidence_directness"].fillna("-"),
        "Evidence basis": basis,
        "Citation": cites,
        "Predicted DDG, isolated subunit (kcal/mol)": lookup(pred, "pred_monomer"),
        "Predicted DDG, in complex (kcal/mol)": lookup(pred, "pred_fold"),
        "Predicted DDG, binding (kcal/mol)": lookup(pred, "pred_binding"),
        "Measured DDG (kcal/mol)": meas,
        "Measured axis": merged["measured_axis"].fillna("-"),
        "Structure": merged["structure_source"].fillna("-"),
    })
    return out.sort_values(["System", "Variant"]).reset_index(drop=True)


def system_table(canonical):
    frame = canonical.copy()
    frame["system_key"] = frame["system"].map(norm_system)
    graded = frame["mech_consistency_t25"].isin(GRADE_MAP)
    pathogenic = frame["phenotype"].astype(str).str.contains("pathogenic", na=False)

    rows = []
    for key, group in frame.groupby("system_key"):
        mask = frame["system_key"] == key
        structures = [s for s in dict.fromkeys(group["structure_source"].dropna())]
        rows.append({
            "System": SYSTEM_LABEL.get(key, key),
            "Variant gene": "/".join(sorted(set(str(g).upper() for g in group["gene"].dropna()))) or "-",
            "Variants": int(len(group)),
            "Mechanism-gradeable": int((mask & graded).sum()),
            # Named "Clinically pathogenic", not "Pathogenic". It is a clinical
            # annotation counted over EVERY variant in the system, on a
            # different axis from the structural ground truth and independent
            # of it -- which is the paper's thesis, not an inconsistency. It can
            # therefore exceed Mechanism-gradeable, as it does for VHL-elongin C
            # (3 of 3 clinically pathogenic, 2 mechanism-gradeable), because a
            # variant can be clinically pathogenic by a mechanism the supplied
            # structure cannot represent. A column headed "Pathogenic" beside a
            # structural count invites exactly the misreading a lab reader had.
            "Clinically pathogenic": int((mask & pathogenic).sum()),
            "Structure": ", ".join(structures) if structures else "-",
            "Model source": ("Experimental" if any(s[:4].isdigit() or s[0].isdigit()
                                                   for s in structures)
                             else "AlphaFold" if structures else "-"),
        })
    return pd.DataFrame(rows).sort_values("System").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the written tables differ from a fresh build")
    ap.add_argument("--print-only", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    canonical, ledger, calib = load()
    variants = variant_table(canonical, ledger, calib)
    systems = system_table(canonical)

    assert len(systems) == EXPECTED_SYSTEMS, \
        "expected %d systems, built %d" % (EXPECTED_SYSTEMS, len(systems))
    assert len(variants) == EXPECTED_VARIANTS, \
        "expected %d variants, built %d" % (EXPECTED_VARIANTS, len(variants))
    assert int(systems["Variants"].sum()) == EXPECTED_VARIANTS, \
        "system variant counts sum to %d" % int(systems["Variants"].sum())

    # Every blank cell in Table S1 must be blank for a stated reason. A reviewer
    # reading a dash cannot tell "not applicable" from "we never looked", so the
    # three reasons are asserted here and stated in the caption. A NEW blank for
    # any other reason fails the build rather than shipping as an unexplained
    # hole.
    def blank(col):
        return variants[col].astype(str).str.strip().isin(("-", "nan", ""))

    single_chain = variants["Gated interface partners"].isna()
    for col in ("Predicted DDG, in complex (kcal/mol)",
                "Predicted DDG, binding (kcal/mol)"):
        assert blank(col).equals(single_chain), (
            "%s is blank on rows other than the single-chain system: %s"
            % (col, variants.loc[blank(col) & ~single_chain, "Variant"].tolist()))

    uncommitted = variants["Committed axes"].eq("none committed")
    assert blank("Evidence class").equals(uncommitted), (
        "Evidence class is blank on rows that DO commit an axis: %s"
        % variants.loc[blank("Evidence class") & ~uncommitted, "Variant"].tolist())

    assert not blank("Evidence basis").any(), (
        "no evidence basis for %s"
        % variants.loc[blank("Evidence basis"), "Variant"].tolist())
    assert not blank("Predicted DDG, isolated subunit (kcal/mol)").any(), \
        "missing isolated-subunit prediction"

    # Measured energies are sparse by nature, but a variant whose evidence class
    # includes E1 (quantitative energetic) and whose measurement is an ENERGY
    # should have one. The exceptions are variants whose E1 evidence is a
    # dissociation constant rather than a folding energy; those carry the K_D in
    # the basis column, and the set is named so it cannot grow silently.
    KD_ONLY = {"I62V", "R53H", "R78G", "A1381T", "R1334Q"}
    e1 = variants["Evidence class"].astype(str).str.contains("E1")
    missing = set(variants.loc[e1 & blank("Measured DDG (kcal/mol)"), "Variant"])
    assert missing == KD_ONLY, (
        "E1 variants without a measured energy changed: %s"
        % sorted(missing.symmetric_difference(KD_ONLY)))

    print("Table 1: %d systems, %d variants, %d mechanism-gradeable"
          % (len(systems), int(systems["Variants"].sum()),
             int(systems["Mechanism-gradeable"].sum())))
    # Anchors come from the 63-row measured-energy cohort, which is WIDER than
    # the 61 scored variants (it carries external SKEMPI points and the
    # haemoglobin W37 substitution series). Only in-benchmark rows join here.
    print("Table S1: %d variants | %d with evidence basis | %d with a measured anchor"
          % (len(variants),
             int((variants["Evidence basis"] != "-").sum()),
             int((variants["Measured DDG (kcal/mol)"] != "-").sum())))
    labile = int(variants["Same call at every threshold"].eq("no").sum())
    print("           %d of %d graded variants change their call somewhere in "
          "the threshold range" % (labile, int(variants["Same call at every "
                                                        "threshold"].ne("-").sum())))

    if args.print_only:
        print()
        print(systems.to_string(index=False))
        return 0

    paths = {"COMAVI_Table1_systems.csv": systems,
             "COMAVI_TableS1_variants.csv": variants}

    if args.check:
        stale = []
        for name, fresh in paths.items():
            path = args.out_dir / name
            if not path.is_file():
                stale.append("%s (missing)" % name)
                continue
            # Compare the CSVs, not a CSV against an in-memory frame. The
            # artifact is the file, and a round trip is lossy: a formatted
            # "0.00" is re-read as the float 0.0, so comparing the written file
            # against the unwritten frame reported drift on a table that
            # round-trips perfectly. Both sides go through read_csv.
            reread = pd.read_csv(io.StringIO(fresh.to_csv(index=False)))
            if not pd.read_csv(path).astype(str).equals(reread.astype(str)):
                stale.append(name)
        if stale:
            print("FAIL: %s differ(s) from a fresh build" % ", ".join(stale))
            return 1
        print("PASS: cohort tables match a fresh build")
        return 0

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in paths.items():
        frame.to_csv(args.out_dir / name, index=False)
        print("wrote %s" % (args.out_dir / name).relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
