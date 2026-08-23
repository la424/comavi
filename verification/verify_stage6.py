#!/usr/bin/env python3
"""
verify_stage6.py — Downstream-only verification for the v5 framework.

Re-runs Stage 6 (mechanism reclassification + structural_agreement +
mech_consistency + diagnostic columns) against the cached intermediate
file, without touching FoldX or the structural metrics computation.

Compares output to expected v5 headline numbers.

Usage:
  python verify_stage6.py --intermediate /path/to/comavi_v7_results_with_nbhd.csv
"""
import argparse
import sys
from pathlib import Path
import subprocess
import pandas as pd

# v5 expected headlines (from methods_metrics_sketch_v3 final values)
EXPECTED_V5 = {
    'structural_agreement': {
        't10': 0.658, 't15': 0.718, 't20': 0.709,
        't25': 0.718, 'tSAP': 0.709,
    },
    'mech_consistency': {
        't10': 0.602, 't15': 0.693, 't20': 0.716,
        't25': 0.761, 'tSAP': 0.750,
    },
    'threshold_stable': 28,
    'level1_TPR_t10_COMAVI_full': 0.913,
    'level1_TPR_t10_monomer_only': 0.391,
    'level3_HBB_pearson_r': 0.89,
    'level2_pathogenic_detection_t10': 24,
    'level2_pathogenic_gof_detection_t10': 5,
    'level2_benign_silent_t25': 7,
}

# Track C — canonical 61-variant release (reference_outputs/scored_61var_canonical.csv).
#
# These are the numbers the manuscript reports, and they are NOT covered by
# EXPECTED_V5 above: that block dates from the 44-variant generation, whose
# per-axis ground truth predates the tightened evidence standard (ledger §2-§4).
# 14 of those 44 variants had ground-truth tokens changed in re-curation, so a
# Track B pass says nothing about the released canonical numbers. Values below
# are post the v7.1 `structurally_uncommitted` grading correction (ledger §15).
# v7.3 RE-FREEZE. The BRCA1-BRCT monomer-fold cohort is now POOLED into both
# headlines (scripts/apply_brct_pooling_v73.py). Graded rows 47 -> 57: the
# cohort contributes 10 (12 rows less R1699L/R1699Q, whose curated mechanism —
# loss of a phospho-peptide binding site — is not observable on the isolated
# tandem-BRCT structure; these are the same 2 the published measured-vs-FoldX
# correlation excludes, via role_in_cohort == 'fold_intact_function_lost').
#
# Structural agreement 131 -> 133 all-row outputs: +12 monomer axes from the
# cohort, less 1 (R1699L, whose monomer CI is indistinguishable from 0). The
# cohort contributes NO tier axis — the tier is undefined without a partner
# chain, and grading it was scoring a NaN as "did not fire" (v7.3 fix (d)).
EXPECTED_CANONICAL = {
    'structural_agreement': {
        't10': (94, 133), 't15': (99, 133), 't20': (99, 133),
        't25': (99, 133), 'tSAP': (97, 133),
    },
    'mech_consistency': {
        't10': 0.579, 't15': 0.649, 't20': 0.684,
        't25': 0.719, 'tSAP': 0.711,
    },
    'mech_graded_n': 57,
    # Four-way decomposition of the t=2.5 headline (ledger §13).
    # monomer 16 -> 27 gradeable: +12 cohort rows, -1 CI-excluded (R1699L).
    'axis_decomposition_t25': {
        'monomer': (21, 28), 'fold': (20, 26),
        'binding': (24, 32), 'tier': (34, 47),
    },
    # Tier pathogenicity gradient — carries no expected_mech_class term, so
    # neither the v7.1 correction nor v7.3 pooling may move it.
    'tier_gradient': [(14, 14), (13, 18), (7, 10), (3, 7)],
    'n_ppi_rows': 49,
    'n_brct_rows': 12,
}

TOLERANCE = 0.005  # within 0.5 percentage points


def check(label, actual, expected, tol=TOLERANCE):
    """Print a checkmark or cross for a numeric assertion."""
    if expected is None:
        print(f"  [SKIP] {label}: no expected value")
        return None
    if actual is None:
        print(f"  [FAIL] {label}: actual is None")
        return False
    delta = abs(actual - expected)
    if delta <= tol:
        print(f"  [ OK ] {label}: {actual:.3f} (expected {expected:.3f})")
        return True
    print(f"  [FAIL] {label}: {actual:.3f} (expected {expected:.3f}, diff {delta:.3f})")
    return False


def run_concordance(input_csv, output_dir, scripts_dir, am_xlsx):
    """Run apply_concordance_v5.py and return path to concordance.csv."""
    cmd = [
        sys.executable,
        str(scripts_dir / 'apply_concordance_v5.py'),
        '--results', str(input_csv),
        '--external', str(am_xlsx),
        '--outdir', str(output_dir),
    ]
    print(f"\n  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\n  ERROR: apply_concordance_v5 failed")
        print(result.stdout[-2000:])
        print("STDERR:")
        print(result.stderr[-2000:])
        sys.exit(1)
    return output_dir / 'comavi_v7_concordance.csv'


def verify_canonical(canonical_csv, scripts_dir):
    """Track C — verify the released canonical 61-variant table.

    There is no 61-variant *pipeline input* in the repo (the canonical table is
    itself the released product), so this track cannot re-run the pipeline end
    to end. Instead it re-derives the grading columns from the released table
    using the pipeline's own functions and checks three things:

      1. the re-derived grades reproduce the STORED columns exactly (guards
         against the table drifting from the code that produced it);
      2. the pooled headlines match EXPECTED_CANONICAL;
      3. the tier gradient is unchanged (it must not depend on grading).

    Two conventions are required and are easy to get wrong (ledger §15):
      - discover_partners() over the wide released table also returns the
        per-partner `_ci95_*` / `_distinguishable_*` columns as if they were
        partner chains; they must be filtered out.
      - grading applies to the 49 interaction rows only; the 12 BRCA1-BRCT
        fold-expansion rows are excluded from interaction-axis grading.
    """
    sys.path.insert(0, str(scripts_dir))
    import apply_concordance_v5 as ac

    print("\n" + "=" * 70)
    print("Track C verification (released canonical 61-variant table)")
    print("=" * 70)

    df = pd.read_csv(canonical_csv, low_memory=False)
    results = []
    MC_MAP = {'consistent': 1.0, 'partial': 0.5, 'inconsistent': 0.0}
    tags = [t for t, _ in ac.THRESHOLD_SPECS]

    partners = [p for p in ac.discover_partners(df)
                if '_ci95_' not in p and '_distinguishable_' not in p]
    # v7.3: partition on the SYSTEM, not on expected_mech_class. The old proxy
    # ("a row is an interaction row iff its expected class is populated") held
    # only while the BRCT cohort's class was unfilled; pooling populated it and
    # silently reclassified all 12 rows as interaction rows.
    brct = df['system'].eq('brca1_brct')
    ppi = ~brct

    print(f"\n  rows: {len(df)}  PPI: {int(ppi.sum())}  BRCT: {int((~ppi).sum())}"
          f"  partner chains: {len(partners)}")
    results.append(check("n_ppi_rows", float(ppi.sum()),
                         float(EXPECTED_CANONICAL['n_ppi_rows']), tol=0))
    results.append(check("n_brct_rows", float((~ppi).sum()),
                         float(EXPECTED_CANONICAL['n_brct_rows']), tol=0))

    # --- 0. structural-provenance completeness (v7.4 regression guard) ------
    #
    # The v7.4 defect: the BRCT supplement was merged into the canonical table
    # WITHOUT site_plddt_status / structure_evaluable, so those 12 rows fell
    # silently out of the concordance framework and the confidence-restriction
    # analysis. apply_concordance_v5 reads site_plddt_status but never writes
    # it, so nothing downstream noticed. Assert completeness on every row so a
    # future supplement merge that drops these columns fails loudly here.
    print("\n  structural-provenance completeness (v7.4 guard):")
    for col, allowed in (("site_plddt_status",
                          {"crystal", "confident", "partial", "low"}),
                         ("structure_evaluable", None)):
        if col not in df.columns:
            print(f"    [FAIL] {col}: column absent from canonical table")
            results.append(False)
            continue
        n_missing = int(df[col].isna().sum())
        results.append(check(f"{col}_complete", float(n_missing), 0.0, tol=0))
        if allowed is not None:
            bad = sorted(set(df[col].dropna().astype(str).str.lower()) - allowed)
            if bad:
                print(f"    [FAIL] {col}: unexpected values {bad}")
                results.append(False)
            else:
                print(f"    [PASS] {col}: all values in {sorted(allowed)}")
                results.append(True)

    # --- 1. re-derive grading from the pipeline's own functions -------------
    #
    # v7.3: the BRCT cohort is now graded too, so the re-derivation must cover
    # ALL rows. Both cohort rules live in the shipped module — the direction
    # -aware fold branch inside derive_expected_mech_class(), and the
    # unobservable-mechanism exclusion in unobservable_variants() — so this
    # track exercises the released code path rather than a local restatement
    # of it. It still re-derives every column independently of the applier.
    unobservable = ac.unobservable_variants()
    reder = df.copy()
    reder['expected_mech_class'] = reder.apply(
        ac.derive_expected_mech_class, axis=1)
    gradeable = ~reder['variant'].isin(unobservable)
    axes_by_idx = {i: ac.classify_axis_status(r) for i, r in reder.iterrows()}
    for suf in tags:
        reder[f'mech_consistency_{suf}'] = None
        grades = [
            ac.grade_mechanism_consistency(
                r, r.get(f'mech_{suf}'), r.get('expected_mech_class'),
                axes_by_idx.get(r.name))[0]
            for _, r in reder[gradeable].iterrows()
        ]
        reder.loc[gradeable, f'mech_consistency_{suf}'] = grades
    for tag, t in ac.THRESHOLD_SPECS:
        mt, ft, bt = ((t['monomer'], t['fold'], t['binding'])
                      if isinstance(t, dict) else (float(t),) * 3)
        sa = reder.apply(
            lambda r, _m=mt, _f=ft, _b=bt:
                ac.compute_structural_agreement(r, partners, _m, _f, _b), axis=1)
        reder[f'structural_agreement_n_{tag}'] = [x[0] for x in sa]
        reder[f'structural_agreement_d_{tag}'] = [x[1] for x in sa]

    print("\n  re-derived grading reproduces stored columns:")
    for tag in tags:
        sn = (int(df[f'structural_agreement_n_{tag}'].fillna(0).sum()),
              int(df[f'structural_agreement_d_{tag}'].fillna(0).sum()))
        rn = (int(reder[f'structural_agreement_n_{tag}'].fillna(0).sum()),
              int(reder[f'structural_agreement_d_{tag}'].fillna(0).sum()))
        ok = sn == rn
        print(f"    [{' OK ' if ok else 'FAIL'}] SA_{tag} stored {sn[0]}/{sn[1]}"
              f" vs re-derived {rn[0]}/{rn[1]}")
        results.append(ok)
        smc = df[f'mech_consistency_{tag}'].map(MC_MAP).dropna()
        rmc = reder[f'mech_consistency_{tag}'].map(MC_MAP).dropna()
        ok = (len(smc) == len(rmc)) and abs(smc.mean() - rmc.mean()) < 1e-9
        print(f"    [{' OK ' if ok else 'FAIL'}] MC_{tag} stored {smc.mean():.4f}"
              f" (n={len(smc)}) vs re-derived {rmc.mean():.4f} (n={len(rmc)})")
        results.append(ok)

    # --- 2. pooled headlines ------------------------------------------------
    print("\n  structural_agreement vs expected:")
    for tag in tags:
        en, ed = EXPECTED_CANONICAL['structural_agreement'][tag]
        n = int(df[f'structural_agreement_n_{tag}'].fillna(0).sum())
        d = int(df[f'structural_agreement_d_{tag}'].fillna(0).sum())
        ok = (n, d) == (en, ed)
        print(f"    [{' OK ' if ok else 'FAIL'}] SA_{tag}: {n}/{d} = {n/d:.4f}"
              f" (expected {en}/{ed} = {en/ed:.4f})")
        results.append(ok)

    print("\n  mech_consistency vs expected:")
    graded = df['mech_consistency_t25'].map(MC_MAP).dropna()
    results.append(check("mech_graded_n", float(len(graded)),
                         float(EXPECTED_CANONICAL['mech_graded_n']), tol=0))
    for tag in tags:
        vals = df[f'mech_consistency_{tag}'].map(MC_MAP).dropna()
        results.append(check(f"mech_consistency_{tag}", float(vals.mean()),
                             EXPECTED_CANONICAL['mech_consistency'][tag]))

    # --- 2b. four-way axis decomposition of the t=2.5 headline --------------
    # v7.3: this constant was previously declared and never asserted. It is the
    # decomposition quoted in the ledger and in the manuscript, so it must be
    # checked, and it must close to the headline it decomposes.
    print("\n  axis decomposition of SA_t25 (must sum to the headline):")
    spec = dict(ac.THRESHOLD_SPECS)['t25']
    triple = ((spec['monomer'], spec['fold'], spec['binding'])
              if isinstance(spec, dict) else (spec, spec, spec))
    dec = {k: [0, 0] for k in ('monomer', 'fold', 'binding', 'tier')}
    for _, r in df.iterrows():
        for k, (a, b) in ac.structural_agreement_by_axis(r, partners, *triple).items():
            dec[k][0] += a
            dec[k][1] += b
    for k in ('monomer', 'fold', 'binding', 'tier'):
        got = tuple(dec[k])
        exp = tuple(EXPECTED_CANONICAL['axis_decomposition_t25'][k])
        ok = got == exp
        print(f"    [{' OK ' if ok else 'FAIL'}] {k}: {got[0]}/{got[1]}"
              f"  (expected {exp[0]}/{exp[1]})")
        results.append(ok)
    dsum = (sum(v[0] for v in dec.values()), sum(v[1] for v in dec.values()))
    hsum = tuple(EXPECTED_CANONICAL['structural_agreement']['t25'])
    ok = dsum == hsum
    print(f"    [{' OK ' if ok else 'FAIL'}] decomposition closes: "
          f"{dsum[0]}/{dsum[1]} vs headline {hsum[0]}/{hsum[1]}")
    results.append(ok)

    # --- 3. tier gradient (must be untouched by grading) --------------------
    print("\n  tier pathogenicity gradient (independent of expected_mech_class):")
    sub = df[df['comavi_tier'].notna() & df['role'].notna()].copy()
    sub['path'] = sub['role'].astype(str).str.startswith('pathogenic')
    for i, tier in enumerate(['Tier 1', 'Tier 2', 'Tier 3', 'Tier 4']):
        rows = sub[sub['comavi_tier'].astype(str) == tier]
        got = (int(rows['path'].sum()), int(len(rows)))
        exp = tuple(EXPECTED_CANONICAL['tier_gradient'][i])
        ok = got == exp
        pct = 100 * got[0] / got[1] if got[1] else 0
        print(f"    [{' OK ' if ok else 'FAIL'}] {tier}: {got[0]}/{got[1]}"
              f" ({pct:.0f}%)  expected {exp[0]}/{exp[1]}")
        results.append(ok)

    return results


def run_evaluation(corrected_csv, output_dir, scripts_dir):
    """Run run_evaluate.py to produce Track A outputs."""
    cmd = [
        sys.executable,
        str(scripts_dir / 'run_evaluate.py'),
        '--input', str(corrected_csv),
        '--output-dir', str(output_dir),
    ]
    print(f"\n  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\n  ERROR: run_evaluate failed")
        print(result.stdout[-2000:])
        print("STDERR:")
        print(result.stderr[-2000:])
        sys.exit(1)
    return output_dir / 'evaluation'


def verify_concordance(concordance_csv):
    """Verify Track B (concordance) headlines."""
    df = pd.read_csv(concordance_csv)
    results = []
    print("\n" + "=" * 70)
    print("Track B verification (apply_concordance_v5)")
    print("=" * 70)

    # structural_agreement
    print("\nstructural_agreement at all 5 thresholds:")
    for tag in ('t10', 't15', 't20', 't25', 'tSAP'):
        n_col = f'structural_agreement_n_{tag}'
        d_col = f'structural_agreement_d_{tag}'
        if n_col not in df.columns or d_col not in df.columns:
            print(f"  [SKIP] {tag}: columns not present")
            continue
        n = df[n_col].fillna(0).sum()
        d = df[d_col].fillna(0).sum()
        score = float(n / d) if d > 0 else 0.0
        ok = check(f"structural_agreement_{tag}", score,
                   EXPECTED_V5['structural_agreement'].get(tag))
        results.append(ok)

    # mech_consistency
    print("\nmech_consistency at all 5 thresholds:")
    for tag in ('t10', 't15', 't20', 't25', 'tSAP'):
        col = f'mech_consistency_{tag}'
        if col not in df.columns:
            print(f"  [SKIP] {tag}: column not present")
            continue
        vc = df[col].value_counts()
        n = len(df) - df[col].isna().sum()
        score = float((vc.get('consistent', 0) + 0.5 * vc.get('partial', 0)) / n) if n else 0.0
        ok = check(f"mech_consistency_{tag}", score,
                   EXPECTED_V5['mech_consistency'].get(tag))
        results.append(ok)

    # threshold-stable
    if 'mech_consistency_threshold_stable' in df.columns:
        n_stable = int(df['mech_consistency_threshold_stable'].sum())
        ok = check("threshold_stable_count", n_stable,
                   EXPECTED_V5['threshold_stable'], tol=2)  # allow ±2 (small variation OK)
        results.append(ok)

    # Spot-check worked example: PIK3CA E545K
    print("\nWorked example: PIK3CA E545K")
    sub = df[(df['gene'].str.lower() == 'pik3ca') & (df['variant'] == 'E545K')]
    if len(sub) > 0:
        r = sub.iloc[0]
        for tag in ('t10', 'tSAP'):
            mech = r[f'mech_{tag}']
            grade = r[f'mech_consistency_{tag}']
            print(f"  {tag}: mech='{mech}', consistency={grade}")

    return results


def verify_evaluation(eval_dir):
    """Verify Track A headlines."""
    print("\n" + "=" * 70)
    print("Track A verification (run_evaluate)")
    print("=" * 70)
    results = []

    # Level 1
    l1_csv = eval_dir / 'level1_merged.csv'
    if l1_csv.exists():
        l1 = pd.read_csv(l1_csv)
        full = l1[l1['classifier'] == 'COMAVI_full']
        mono = l1[l1['classifier'] == 'monomer_only']
        if len(full) > 0:
            ok = check("Level 1 COMAVI_full TPR @ t=1.0",
                       float(full.iloc[0]['TPR_t10']),
                       EXPECTED_V5['level1_TPR_t10_COMAVI_full'], tol=0.005)
            results.append(ok)
        if len(mono) > 0:
            ok = check("Level 1 monomer_only TPR @ t=1.0",
                       float(mono.iloc[0]['TPR_t10']),
                       EXPECTED_V5['level1_TPR_t10_monomer_only'], tol=0.005)
            results.append(ok)

    # Level 3 HBB
    l3_csv = eval_dir / 'level3_hbb_quantitative.csv'
    if l3_csv.exists():
        l3 = pd.read_csv(l3_csv)
        if len(l3) > 0:
            ok = check("Level 3 HBB Pearson r",
                       float(l3.iloc[0]['pearson_r']),
                       EXPECTED_V5['level3_HBB_pearson_r'], tol=0.01)
            results.append(ok)

    # Level 2 per-phenotype detection
    l2_csv = eval_dir / 'level2_merged_detection.csv'
    if l2_csv.exists():
        l2 = pd.read_csv(l2_csv)
        for phen, key in [
            ('pathogenic', 'level2_pathogenic_detection_t10'),
            ('pathogenic_gof', 'level2_pathogenic_gof_detection_t10'),
        ]:
            sub = l2[l2['phenotype'] == phen]
            if len(sub) > 0:
                actual = int(sub.iloc[0]['n_detected_t10'])
                ok = check(f"Level 2 {phen} detected @ t=1.0",
                           actual, EXPECTED_V5[key], tol=1)
                results.append(ok)

    return results


def main():
    parser = argparse.ArgumentParser(description="Verify v5 framework downstream")
    parser.add_argument('--intermediate', required=True,
                        help='Path to comavi_v7_results_with_nbhd.csv')
    parser.add_argument('--am', required=True,
                        help='Path to AM xlsx')
    parser.add_argument('--scripts-dir', required=True,
                        help='Path to PIPELINE_CURRENT/scripts/')
    parser.add_argument('--corrected',
                        help='Path to comavi_v7_results_corrected.csv (for Track A)')
    parser.add_argument('--canonical',
                        help='Path to reference_outputs/scored_61var_canonical.csv '
                             '(Track C — the released manuscript numbers)')
    parser.add_argument('--output-dir', default='./verification_output',
                        help='Where to write verification artifacts')
    args = parser.parse_args()

    intermediate = Path(args.intermediate).resolve()
    am_xlsx = Path(args.am).resolve()
    scripts_dir = Path(args.scripts_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(exist_ok=True, parents=True)

    print(f"Verifying v5 framework using:")
    print(f"  intermediate: {intermediate}")
    print(f"  AM xlsx:      {am_xlsx}")
    print(f"  scripts:      {scripts_dir}")
    print(f"  output:       {output_dir}")

    # Track B
    concordance_csv = run_concordance(intermediate, output_dir, scripts_dir, am_xlsx)
    track_b_results = verify_concordance(concordance_csv)

    # Track A — needs the corrected CSV (output of baseline_correct, not concordance)
    track_a_results = []
    if args.corrected:
        eval_dir = run_evaluation(Path(args.corrected).resolve(), output_dir, scripts_dir)
        track_a_results = verify_evaluation(eval_dir)
    else:
        print("\n[SKIP Track A] --corrected not provided")

    # Track C — released canonical 61-variant table (the manuscript numbers)
    track_c_results = []
    if args.canonical:
        track_c_results = verify_canonical(Path(args.canonical).resolve(), scripts_dir)
    else:
        print("\n[SKIP Track C] --canonical not provided")

    # Summary
    all_results = [r for r in (track_b_results + track_a_results + track_c_results)
                   if r is not None]
    n_pass = sum(1 for r in all_results if r is True)
    n_fail = sum(1 for r in all_results if r is False)
    print("\n" + "=" * 70)
    print(f"Verification summary: {n_pass}/{n_pass + n_fail} passed")
    print("=" * 70)
    if n_fail > 0:
        print(f"\n[FAIL] {n_fail} checks failed. Review output above.")
        sys.exit(1)
    print("\n[PASS] All v5 framework checks passed.")
    sys.exit(0)


if __name__ == '__main__':
    main()
