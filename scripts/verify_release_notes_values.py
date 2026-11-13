#!/usr/bin/env python3
"""Assert the current-column values in RELEASE_NOTES_*.md match the records.

WHY THIS EXISTS
---------------
Release notes must restate headline numbers -- that is what they are for -- so
verify_no_duplicated_headlines.py cannot simply forbid them. But an unchecked
document that restates a headline is exactly how this project has produced drift
before: a figures README quoting a cohort that no longer existed, a results
summary whose template hardcoded a fraction beside the computed rate.

So the notes are allowed there and gated here instead. Every needle below is
BUILT from a committed record and asserted to appear in the notes, in the bold
form the current column uses. If a record moves and the notes are not updated,
this fails. The previous-release column is deliberately not checked against the
records -- those values are superseded by construction and are verified against
what the previous tag actually published, which git holds and these notes
document.

v11: a release's notes describe that release. Until v11 every notes file was
held to the working-tree records, which only worked while no value moved after
a release: the first value-changing revision after v3.0-methods made the v3.0
notes "fail" for stating exactly what v3.0 published. A notes file whose
version is a git tag (RELEASE_NOTES_<tag>.md) is now held to the records AT
THAT TAG, read with git show; a notes file without a tag is a release in
preparation (as verify_citation_version.py already allows) and is held to the
working tree. Values that moved since the newest tag are listed, so the notes
for the next release know what they must restate. --require-history makes a
missing tag a failure; CI passes it, with a full-history checkout.

Usage:
    python scripts/verify_release_notes_values.py [--require-history]
    python scripts/verify_release_notes_values.py --notes RELEASE_NOTES_v3.0-methods.md
"""
import argparse
import glob
import io
import json
import pathlib
import re
import subprocess
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
RO = REPO / "reference_outputs"


def reader(rev=None):
    """Read a repository file from the working tree (rev None) or a git revision."""
    if rev is None:
        return lambda rel: (REPO / rel).read_text()

    def read(rel):
        out = subprocess.run(["git", "-C", str(REPO), "show", "%s:%s" % (rev, rel)],
                             capture_output=True, text=True)
        if out.returncode:
            raise FileNotFoundError("%s:%s -- %s" % (rev, rel, out.stderr.strip()[:120]))
        return out.stdout
    return read


def tag_exists(tag):
    return subprocess.run(["git", "-C", str(REPO), "rev-parse", "-q", "--verify", "refs/tags/" + tag],
                          capture_output=True).returncode == 0


def expected_current(read=None):
    """[(label, needle)] -- every needle computed, none typed."""
    read = read or reader()
    sweep = {r["threshold_tag"]: r for r in
             json.loads(read("reference_outputs/COMAVI_threshold_sweep.json"))["rows"]}
    t25 = sweep["t25"]
    strat = json.loads(read("reference_outputs/COMAVI_evidence_stratified_agreement.json"))
    stress = json.loads(read("reference_outputs/stress_tests/COMAVI_stress_verified_statistics.json"))
    led = pd.read_csv(io.StringIO(read("reference_outputs/COMAVI_evidence_ledger.csv")))
    pm = pd.read_csv(io.StringIO(read("reference_outputs/isds_v1/ISDS_v1_primary_metrics.csv"))).set_index("score")
    pv = pd.read_csv(io.StringIO(read("reference_outputs/isds_v1/ISDS_v1_per_variant.csv")))
    can = pd.read_csv(io.StringIO(read("reference_outputs/scored_61var_canonical.csv")), low_memory=False)
    ranked = (pv.sort_values("isds_v1", ascending=False)
                ["structural_ground_truth"].astype(bool).values)
    npos = int(ranked.sum())
    ci = stress["system_cluster_ci95"]["mc"]["reportable_3"]
    c, acc, st = strat["crude"], strat["accounting"], strat["stratified"]
    e1 = strat["by_evidence_type"]["E1_quantitative_energetic"]

    # (label, needle, in_comparison_table). Table rows must appear in the bold
    # current column; a value present only unbolded there is the superseded one.
    return [
        ("mechanism-pattern agreement",
         "%.4f (%.1f/%d)" % (t25["whole_variant"], t25["whole_variant_k"],
                             t25["whole_variant_n"])),
        ("graded for whole-variant localization",
         "%d of %d" % (t25["whole_variant_n"], len(can))),
        ("committed axis expectations", "%d rows" % len(led)),
        ("four-output continuity aggregate",
         "%d/%d = %.3f" % (t25["structural_agreement_k"],
                           t25["structural_agreement_n"],
                           t25["structural_agreement"])),
        ("prioritization population",
         "%d (%d / %d)" % (len(pv), npos, len(pv) - npos)),
        ("priority-score ROC AUC", "%.3f" % pm.loc["isds_v1", "roc_auc"]),
        ("top-20 recovery", "%d/%d" % (int(ranked[:20].sum()), npos)),
        # Stated in prose rather than in the comparison table, so presence is
        # sufficient: there is no superseded twin on the same row for a bold
        # marker to distinguish it from.
        ("system-cluster 95% interval", "%.3f\u2013%.3f" % (ci[0], ci[1]), False),
        # Supplementary evidence stratification
        ("committed axes / total gradeable",
         "%d / %d" % (acc["n_committed_axes"], acc["n_total_gradeable"])),
        ("E1\u2013E2 quantitative evidence",
         "%d/%d = %.3f" % (c["quant_ok"], c["quant_n"], c["quant_rate"])),
        ("E3\u2013E5 qualitative or inferred",
         "%d/%d = %.3f" % (c["soft_ok"], c["soft_n"], c["soft_rate"])),
        ("crude contrast, Fisher p", "%.4f" % c["p"]),
        ("Mantel-Haenszel OR",
         "%.3f (p = %.3f)" % (st["OR_mantel_haenszel"], st["p"])),
        ("E1 quantitative energetic axes",
         "%d/%d = %.3f" % (e1["ok"], e1["n"], e1["rate"])),
    ]


def notes_version(path):
    m = re.match(r"RELEASE_NOTES_(.+)\.md$", path.name)
    return m.group(1) if m else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--notes", default=None,
                    help="a single notes file; default is every RELEASE_NOTES_*.md")
    ap.add_argument("--require-history", action="store_true",
                    help="fail, rather than skip, when a released version's tag is not in the clone")
    args = ap.parse_args()
    files = ([REPO / args.notes] if args.notes
             else sorted(pathlib.Path(p) for p in glob.glob(str(REPO / "RELEASE_NOTES_*.md"))))
    if not files:
        print("FAIL: no RELEASE_NOTES_*.md found; this gate is wired into CI and "
              "a zero-file glob must not report success")
        return 1

    bad, covered, checked, skipped = [], set(), [], []
    needles_by_rev = {}
    for path in files:
        version = notes_version(path)
        rev = version if version and tag_exists(version) else None
        if rev is None and version and _older_than_citation(version):
            # a released version whose tag is not in this clone (shallow checkout
            # or an archive without history)
            if args.require_history:
                bad.append((path.name, "history", version, "tag not found; check out with full history"))
            else:
                skipped.append(path.name)
            continue
        if rev not in needles_by_rev:
            needles_by_rev[rev] = expected_current(reader(rev))
        needles = needles_by_rev[rev]
        checked.append("%s (against %s)" % (path.name, rev or "the working tree"))
        text = path.read_text()
        for label, needle, *rest in needles:
            in_table = rest[0] if rest else True
            # A patch release's notes restate only what that release changed, so
            # a file is held to the rows it actually claims. A CLAIM is the
            # quantity named in a table row, not merely mentioned in prose: the
            # v3.0.1 notes discuss "the mechanism-pattern agreement" without
            # restating it, which is not a claim about its value and must not be
            # treated as one. On its own this would let a row vanish silently, so
            # after the loop every label must be claimed by some notes file.
            if in_table:
                rows = [ln for ln in text.splitlines()
                        if ln.lstrip().startswith("|") and label in ln]
                if not rows:
                    continue
                claim = "\n".join(rows)
            else:
                if label not in text:
                    continue
                claim = text
            covered.add(label)
            # The current column is bold; a value present only unbolded is in the
            # superseded column and does not satisfy the claim.
            if in_table:
                ok = ("**%s**" % needle) in claim or ("**%s" % needle) in claim
                why = ("present but not in the bold current column"
                       if needle in claim else "absent from its table row")
            else:
                ok = needle in claim
                why = "absent"
            if not ok:
                bad.append((path.name, label, needle, why))
    for name, label, needle, why in bad:
        print("  MISSING [%s] %s: expected %r, %s" % (name, label, needle, why))
    labels = [lab for lab, *_ in expected_current()]
    uncovered = [lab for lab in labels if lab not in covered] if not skipped else []
    for lab in uncovered:
        print("  UNCLAIMED [%s]: no release-notes file states this quantity, so "
              "nothing gates it" % lab)
    for name in skipped:
        print("  SKIP [%s]: its release tag is not in this clone, so the values it published cannot be "
              "read; CI checks it with --require-history" % name)
    # what the next release's notes will have to restate
    tagged = [notes_version(p) for p in files if notes_version(p) and tag_exists(notes_version(p))]
    if tagged:
        newest = sorted(tagged, key=_version_key)[-1]
        then = dict((lab, nd) for lab, nd, *_ in needles_by_rev.get(newest) or expected_current(reader(newest)))
        now = dict((lab, nd) for lab, nd, *_ in expected_current())
        moved = [lab for lab in now if then.get(lab) != now[lab]]
        if moved:
            print("  NOTE: %d tracked quantit%s moved since %s; the next release's notes must restate %s: %s"
                  % (len(moved), "y has" if len(moved) == 1 else "ies have", newest,
                     "it" if len(moved) == 1 else "them", "; ".join("%s %s -> %s" % (m, then.get(m), now[m]) for m in moved)))
    if bad or uncovered:
        print("FAIL: %d value(s) disagree with the records of their release, %d quantity(ies) "
              "unclaimed across %d notes file(s)" % (len(bad), len(uncovered), len(files)))
        return 1
    print("PASS: %d release-notes file(s) checked, each against the records of its own release (%s); "
          "every quantity is claimed" % (len(checked), "; ".join(checked)))
    return 0


def _version_key(v):
    """Numeric release order: v3.0-methods < v3.0.1-methods < v3.1-methods."""
    return tuple(int(x) for x in re.findall(r"\d+", v))


def _older_than_citation(version):
    """True when a notes file without a tag in this clone names a release older
    than the one CITATION.cff declares, so its values belong to a tag this clone
    lacks. The notes for the declared version itself are checked against the
    working tree, which is that release in its archive and the release in
    preparation in the repository."""
    cff = (REPO / "CITATION.cff").read_text() if (REPO / "CITATION.cff").exists() else ""
    m = re.search(r'^version:\s*"?([^"\n]+)"?', cff, re.M)
    if not m:
        return False
    declared = m.group(1).strip()
    declared = declared if declared.startswith("v") else "v" + declared
    return _version_key(version) < _version_key(declared)


if __name__ == "__main__":
    sys.exit(main())
