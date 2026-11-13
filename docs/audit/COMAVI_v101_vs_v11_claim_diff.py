"""Claim-by-claim diff of the manuscript text, v10.1 (HEAD) -> v11 (working tree).

Both inputs are the text extracted from the built .docx (docs/COMAVI_manuscript_current.txt),
so the comparison is of what a reader sees. Citation numbers are stripped before comparing,
because v11 renumbers references and a renumbering is not a change of claim.

Sentences are aligned in document order on their wording with every number masked, so a
sentence whose wording survives but whose numbers moved is reported as a number change,
not as a deletion plus an insertion. Sentences left unpaired in a changed block are then
matched across the whole document, which catches passages that moved section.
"""
import collections
import difflib
import re

CIT = re.compile(r"\s?\[\d+(?:\s*[,\u2013\-]\s*\d+)*\]")
NUM = re.compile(r"(?<![A-Za-z0-9.])[\u2212\-]?\d+(?:[.,]\d+)*(?:\u00d710\^?[\u2212\-]?\d+)?")
# references and names that contain digits but state no value
REFS = re.compile(r"\b(?:Figs?|Table|Tables)\s+\d+[A-Z]?(?:\s*(?:and|,|\u2013|-)\s*\d*[A-Z])*"
                  r"|\bS\d\s+(?:Text|Table|Fig)\b|\b[A-Z]\d+[A-Z]\b|\b[A-Z]{2,}\d+[A-Z0-9]*\b|\b\d+[A-Z]{2,}\d*\b"
                  r"|\bPDB\s+\w+|\bE[1-5]\b|\bTier\s+\d|\bTiers\s+\d\s*(?:\u2013|-|to)\s*\d|\bAlphaFold\s+3|\bFoldX\s+5\.\d|\bISDS-v1|\bv\d+(?:\.\d+)?")
ABBR = ["e.g.", "i.e.", "et al.", "vs.", "Fig.", "Figs.", "approx.", "cf.", "No.", "Ref."]


def norm(s):
    s = CIT.sub("", s)
    s = s.replace("\u2009", " ").replace("\u00a0", " ")
    return re.sub(r"\s+", " ", s).strip()


def sentences(par):
    p = par
    for i, a in enumerate(ABBR):
        p = p.replace(a, a.replace(".", "\u0000%d" % i))
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(\u0394\u03c1\[])", p)
    out = []
    for s in parts:
        for i, a in enumerate(ABBR):
            s = s.replace(a.replace(".", "\u0000%d" % i), a)
        s = s.strip()
        if s:
            out.append(s)
    return out


def is_heading(line):
    t = line.strip()
    return (0 < len(t) < 130 and not t.endswith((".", ":", ";", ","))
            and "\t" not in t and not re.match(r"^\[|^\d+\.\s", t) and not t.startswith(("Short title", "*")))


def units(text):
    """(part, section, sentence) for the main text, the main tables and the S1-S3 Text."""
    lines = text.split("\n")
    out, part, section, in_refs, in_table = [], "main", "Title", False, None
    for line in lines:
        t = line.strip()
        if not t:
            continue
        if t.startswith("[SI ") and t.endswith("_Table]"):
            break                                     # workbook sheets are checked cell by cell elsewhere
        if t.startswith("[SI "):
            part, section, in_refs, in_table = t[4:-1].replace("_", " "), t[4:-1].replace("_", " "), False, None
            continue
        if t.startswith("[TABLE"):
            in_table = t[1:-1].title()
            continue
        if t == "References":
            in_refs = True
            continue
        if t == "Supporting information":
            in_refs, section = False, "Supporting information"
            continue
        if in_refs:
            continue
        if in_table and "\t" in line:
            out.append((part, "%s, %s" % (section if part != "main" else "Main", in_table),
                        norm(" | ".join(c.strip() for c in line.split("\t")))))
            continue
        in_table = None
        if is_heading(t):
            section = t
            continue
        for s in sentences(norm(t)):
            out.append((part, section, s))
    return out


def skel(s):
    return NUM.sub("#", s)


def nums(s):
    return [m.group(0).replace("\u2212", "-") for m in NUM.finditer(REFS.sub(" ", s))]


def ratio(a, b):
    return difflib.SequenceMatcher(None, a.split(), b.split(), autojunk=False).ratio()


def diff(old_text, new_text, pair_min=0.5, move_min=0.75):
    A, B = units(old_text), units(new_text)
    sa, sb = [skel(x[2]) for x in A], [skel(x[2]) for x in B]
    sm = difflib.SequenceMatcher(None, sa, sb, autojunk=False)
    rows, loose_a, loose_b = [], [], []
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == "equal":
            for i, j in zip(range(i1, i2), range(j1, j2)):
                if A[i][2] != B[j][2]:
                    rows.append(("numbers changed", i, j))
            continue
        ia, jb = list(range(i1, i2)), list(range(j1, j2))
        used = set()
        for i in ia:
            best = max(((ratio(sa[i], sb[j]), j) for j in jb if j not in used), default=(0, None))
            if best[1] is not None and best[0] >= pair_min:
                used.add(best[1])
                rows.append(("pair", i, best[1]))
            else:
                loose_a.append(i)
        loose_b += [j for j in jb if j not in used]
    # passages that moved: match the leftovers across the whole document
    moved_b = set()
    still_a = []
    for i in loose_a:
        best = max(((ratio(sa[i], sb[j]), j) for j in loose_b if j not in moved_b), default=(0, None))
        if best[1] is not None and best[0] >= move_min:
            moved_b.add(best[1])
            rows.append(("moved", i, best[1]))
        else:
            still_a.append(i)
    rows += [("removed", i, None) for i in still_a]
    rows += [("added", None, j) for j in loose_b if j not in moved_b]
    out = []
    for kind, i, j in rows:
        a = A[i] if i is not None else None
        b = B[j] if j is not None else None
        na, nb = collections.Counter(nums(a[2])) if a else collections.Counter(), \
            collections.Counter(nums(b[2])) if b else collections.Counter()
        if kind == "pair":
            kind = "reworded, numbers changed" if na != nb else "reworded"
        if kind == "moved":
            kind = "moved, numbers changed" if na != nb else ("moved, reworded" if a[2] != b[2] else "moved")
        out.append({
            "order_v11": j if j is not None else "", "order_v101": i if i is not None else "",
            "part": (b or a)[0], "section_v101": a[1] if a else "", "section_v11": b[1] if b else "",
            "change": kind, "text_v101": a[2] if a else "", "text_v11": b[2] if b else "",
            "numbers_removed": "; ".join(sorted((na - nb).elements())),
            "numbers_added": "; ".join(sorted((nb - na).elements())),
        })
    key = lambda r: (r["order_v11"] if r["order_v11"] != "" else 10 ** 6, r["order_v101"] if r["order_v101"] != "" else 0)
    return sorted(out, key=key), len(A), len(B)


# ---------------------------------------------------------------- causes
# Families of change between v10.1 and v11. Assigned by the rules below, in order;
# the rule that fired is written beside each row (cause_basis). Rows the rules
# leave unexplained are set by hand in MANUAL, keyed by the start of the v11 text.
VERIFY_MARKS = [  # bounded patterns for the sentences corrected while verifying v11 (draft-stage fixes and
    # docs/audit/COMAVI_v11_post_draft_text_edits.csv); each matches only its own corrected sentence
    r"(?<![0-9])8 at the prediction-interval limit", r"mean of 0\.362", r"replicate mean \u00b1 1\.96",
    r"79 of the 97 expectations qualify", r"Ten binding expectations rest only", r"expectations of a change measured only",
    r"six interface expectations of a change", r"0\.642 to 0\.686", r"structure or conformational ensemble",
    r"scarce and heterogeneous enough to complicate reliable benchmarking",
    r"lowering the binding threshold anywhere in the prespecified range",
    r"Each of the 5 interface mechanisms not localized at the reference",
    r"Ten committed destabilizing expectations fall below even 1\.0"]
KW = [(f, re.compile(p, re.I)) for f, p in (
    ("LABELS", r"AlphaMissense|pathogen|benign|ClinVar|Franklin|clinical label|VUS"),
    ("EVRES", r"amount of complex|either assembled|co-immunoprecipitation|two-hybrid|pull-down|jointly|attribut|localiz|"
              r"axis the literature implicates|wrong axis|right axis"),
    ("RUBRIC", r"\bpartial|inconsistent|\bconsistent|mechanism-pattern score|\bgrade|structure-blind|constant|half credit|\barm\b"),
    ("STATS", r"interval|p = |permut|resampl|bootstrap|Spearman|Fisher|Mann|ROC AUC|average precision"))]
BENCH_OLD, BENCH_NEW = {"61", "99", "49", "12", "47", "57", "98"}, {"56", "97", "46", "10"}
EVRES_NUM = [("77", "79"), ("86", "89"), ("123", "125"), ("15", "18"), ("56", "59"), ("29", "25")]
T1_NOTE = ("Table 1 restructured: the second count column is now the number of structural mechanisms (v10.1: graded "
           "variants); variant counts follow the 56-variant benchmark and the pathogenic column the expert-panel-then-"
           "Franklin label rule")
MANUAL = {
    "(A) The 56 benchmark variants and the populations": ("NEW", "manual: legend of the new cohort-flow panel"),
    "At the 2.5 kcal/mol reference, the isolated-subunit axis alone": (
        "NEW", "manual: the detection gain is now stated with its cost in false calls (27 non-disruptive variants)"),
    "Among the 11 interface variants curated as destabilizers": (
        "FRAMING", "manual: the 12-versus-11 bridge simplified; counts unchanged"),
}
FAMILY = {
    "BENCH": "benchmark restricted to the 56 graded variants",
    "RUBRIC": "grader implements Table 3 as published",
    "EVRES": "evidence resolution for complex-amount expectations",
    "LABELS": "clinical labels: ClinVar expert panel, then Franklin",
    "STATS": "statistic recomputed by a committed generator",
    "VERIFY": "correction found while verifying v11",
    "NEW": "new v11 content",
    "REMOVED": "removed or superseded",
    "FRAMING": "wording only",
    "TERMS": "terminology only",
    "TABLE_RESTRUCTURED": "table restructured",
    "SI_REBUILT": "supplementary text regenerated from the records",
    "RECOMPUTED": "value regenerated, family not identified",
}


def classify(r):
    txt = r["text_v101"] + " " + r["text_v11"]
    rem = set(filter(None, r["numbers_removed"].split("; ")))
    add = set(filter(None, r["numbers_added"].split("; ")))
    kws = [f for f, pat in KW if pat.search(txt)]
    topic = ("; topic: " + ", ".join(kws)) if kws else ""
    if r["part"] != "main":
        return "SI_REBUILT", "S1-S3 Text regenerated from the records by scripts/write_si_text_sources_v11.py" + topic
    for k, v in MANUAL.items():
        if r["text_v11"].startswith(k):
            return v
    hit = [m for m in VERIFY_MARKS if re.search(m, r["text_v11"])]
    if hit:
        return "VERIFY", "v11 verification correction: '%s'" % hit[0]
    if "Table 1" in (r["section_v11"], r["section_v101"]) or r["section_v11"] == "Main, Table 1" \
            or r["section_v101"] == "Main, Table 1":
        return "TABLE_RESTRUCTURED", T1_NOTE
    if "lesion" in r["text_v101"].lower() and "disruption" in r["text_v11"].lower() and not (rem or add):
        return "TERMS", "lesion renamed disruption"
    if r["change"] in ("reworded", "moved, reworded", "moved"):
        return "FRAMING", "wording changed, no value changed"
    if r["change"] == "added":
        return "NEW", "added in v11" + topic
    if r["change"] == "removed":
        return "REMOVED", "removed or superseded in v11" + topic
    if (rem & BENCH_OLD) and (add & BENCH_NEW):
        return "BENCH", "numbers %s -> %s" % ("/".join(sorted(rem & BENCH_OLD)), "/".join(sorted(add & BENCH_NEW)))
    if any(a in rem and b in add for a, b in EVRES_NUM) and re.search(r"attribut|agreement|axis|axes|localiz", txt, re.I):
        return "EVRES", "numbers %s" % ", ".join("%s -> %s" % (a, b) for a, b in EVRES_NUM if a in rem and b in add)
    if kws:
        return kws[0], "keyword: " + ", ".join(kws)
    return "RECOMPUTED", "value regenerated from the records; no family keyword"


def main():
    import argparse
    import csv
    import pathlib
    import subprocess
    ap = argparse.ArgumentParser(description="claim-by-claim diff of the extracted manuscript text")
    ap.add_argument("--repo", required=True)
    ap.add_argument("--base", default="e1dcbb3", help="commit holding the v10.1 text")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    repo = pathlib.Path(a.repo)
    old = subprocess.run(["git", "show", "%s:docs/COMAVI_manuscript_current.txt" % a.base], cwd=repo,
                         capture_output=True, text=True, check=True).stdout
    new = (repo / "docs" / "COMAVI_manuscript_current.txt").read_text()
    rows, n_old, n_new = diff(old, new)
    cols = ["part", "section_v101", "section_v11", "change", "cause", "cause_family", "cause_basis", "text_v101", "text_v11",
            "numbers_removed", "numbers_added", "order_v101", "order_v11"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            c, b = classify(r)
            r.update(cause=c, cause_family=FAMILY[c], cause_basis=b)
            w.writerow({k: r[k] for k in cols})
    print("%d changed units of %d (v10.1) and %d (v11) written to %s" % (len(rows), n_old, n_new, a.out))


if __name__ == "__main__":
    main()
