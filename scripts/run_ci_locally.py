#!/usr/bin/env python3
"""Run locally exactly what the GitHub workflows run, and nothing else.

WHY THIS EXISTS
---------------
Remote CI was red for 20 consecutive runs (2026-10-05 .. 2026-10-21) while local
verification was reported as passing. Both reports were accurate about what they
measured: the local runs executed a hand-picked list of checks, and the failing
check was not on it. A hand-picked list cannot be trusted to represent the
workflow, because the two drift apart silently and in the direction that looks
healthy.

This extracts every `run:` block from every job in every workflow and executes
each block as one shell script, which is what the GitHub runner does. Executing
whole blocks rather than individual lines matters: several blocks contain
multi-line commands with trailing backslashes, and splitting on newlines turns
those into separate, broken invocations.

Steps whose only purpose is provisioning the runner (pip install) are skipped and
reported as skipped rather than silently dropped.

Usage:
    python scripts/run_ci_locally.py                 # every workflow
    python scripts/run_ci_locally.py --workflow verify.yml
    python scripts/run_ci_locally.py --list          # show steps without running
"""
import argparse
import pathlib
import re
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
WORKFLOWS = REPO / ".github" / "workflows"
SKIP = re.compile(r"python -m pip install|pip install --upgrade|apt-get")


def steps(path: pathlib.Path):
    """[(step_name, script)] for every `run:` block, in file order.

    Parsed without PyYAML so this has no dependency the workflows themselves do
    not already have. A `run: |` block is the indented body that follows; a
    single-line `run: cmd` is its own block.
    """
    out, lines = [], path.read_text().splitlines()
    name = "(unnamed)"
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r"\s*-?\s*name:\s*(.+?)\s*$", line)
        if m:
            name = m.group(1).strip("\"'")
            i += 1
            continue
        m = re.match(r"(\s*)-?\s*run:\s*(\|.*|>.*|.*)$", line)
        if m:
            indent, rest = len(m.group(1)), m.group(2).strip()
            if rest.startswith(("|", ">")):
                body, i = [], i + 1
                while i < len(lines):
                    nxt = lines[i]
                    if nxt.strip() and (len(nxt) - len(nxt.lstrip())) <= indent:
                        break
                    body.append(nxt)
                    i += 1
                out.append((name, "\n".join(body)))
            else:
                out.append((name, rest))
                i += 1
            continue
        i += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workflow", action="append",
                    help="file name under .github/workflows (default: all)")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    files = ([WORKFLOWS / w for w in args.workflow] if args.workflow
             else sorted(WORKFLOWS.glob("*.yml")))
    env_line = 'export PYTHONPATH="%s:%s/scripts" MPLBACKEND=Agg\n' % (REPO, REPO)

    passed = failed = skipped = 0
    failures = []
    for wf in files:
        print("=== %s" % wf.name)
        for name, script in steps(wf):
            if not script.strip():
                continue
            if SKIP.search(script):
                print("  SKIP  %-54s (runner provisioning)" % name[:54])
                skipped += 1
                continue
            if args.list:
                first = next((l for l in script.splitlines() if l.strip()), "")
                print("  %-56s %s" % (name[:56], first.strip()[:60]))
                continue
            t0 = time.time()
            proc = subprocess.run(["bash", "-e", "-c", env_line + script],
                                  cwd=REPO, capture_output=True, text=True)
            dt = time.time() - t0
            if proc.returncode == 0:
                print("  PASS  %-54s %5.1fs" % (name[:54], dt))
                passed += 1
            else:
                print("  FAIL  %-54s %5.1fs  rc=%d" % (name[:54], dt, proc.returncode))
                tail = (proc.stdout + proc.stderr).strip().splitlines()
                for line in [l for l in tail if re.search(r"FAIL|Error|error:|Traceback|Assertion", l)][:4]:
                    print("          %s" % line[:150])
                if not any(re.search(r"FAIL|Error", l) for l in tail):
                    for line in tail[-3:]:
                        print("          %s" % line[:150])
                failed += 1
                failures.append("%s / %s" % (wf.name, name))

    if args.list:
        return 0
    print("\n  PASSED %d  FAILED %d  SKIPPED %d" % (passed, failed, skipped))
    for f in failures:
        print("    failing: %s" % f)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
