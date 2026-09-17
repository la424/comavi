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
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
CANON = REPO / "reference_outputs" / "scored_61var_canonical.csv"
LEDGER = REPO / "reference_outputs" / "COMAVI_evidence_ledger.csv"
CALIB = REPO / "reference_outputs" / "COMAVI_delta_calibration_points.csv"
OUT_DIR = REPO / "reference_outputs" / "cohort_tables"

EXPECTED_SYSTEMS = 14
EXPECTED_VARIANTS = 61

GRADE_MAP = {"consistent": 1.0, "partial": 0.5, "inconsistent": 0.0}

# Display names. The canonical's snake_case keys are pipeline identifiers, not
# something a reader should have to decode.
SYSTEM_LABEL = {
    "brca1_bard1": "BRCA1-BARD1 RING heterodimer",
    "brca1_brct": "BRCA1 BRCT tandem domain (fold only)",
    "cam_cav12": "Calmodulin-Ca(v)1.2 IQ domain",
    "cfh_c3b": "Complement factor H-C3b",
    "hemoglobin_dimer": "Haemoglobin alpha-beta dimer",
    "hemoglobin_tetramer": "Haemoglobin tetramer",
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
    "fold_mechanism": "Fold destabilisation",
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

    out = pd.DataFrame({
        "System": merged["system_key"].map(lambda k: SYSTEM_LABEL.get(k, k)),
        "Gene": merged["gene"].fillna("-"),
        "Variant": merged["variant"],
        "Clinical annotation": merged["phenotype"].map(lambda v: pretty(v, PHENO_LABEL)),
        "Expected mechanism": merged["expected_mech_class"].map(lambda v: pretty(v, MECH_LABEL)),
        "Structural tier": merged["comavi_tier"].fillna("-"),
        "COMAVI agreement (t=2.5)": merged["mech_consistency_t25"].fillna("ungraded"),
        "Grantham": merged["grantham_distance"],
        "Gated interface partners": merged["n_interface_partners_gated"],
        "Committed axes": merged["committed_axes"].fillna("none committed"),
        "Evidence class": merged["evidence_class"].fillna("-"),
        "Evidence directness": merged["evidence_directness"].fillna("-"),
        "Evidence basis": merged["evidence_basis"].fillna("-"),
        "Citation": merged["citation"].fillna("-"),
        "Measured DDG (kcal/mol)": merged["measured_kcal"].fillna("-"),
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
            "Pathogenic": int((mask & pathogenic).sum()),
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
            if not pd.read_csv(path).astype(str).equals(fresh.astype(str)):
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
