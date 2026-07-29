#!/usr/bin/env python3
"""
COMAVI brca1_brct monomer-fold FoldX 5.1 BuildModel driver (standalone).
Mirrors scripts/comavi_v7/foldx_runner.build_model: RepairPDB (once) then BuildModel n_runs=5.
Outputs ddg_monomer mean+SD per variant to brct_foldx_ddg.csv.

USAGE:
  python run_brct_foldx.py --foldx /path/to/foldx --rotabase /path/to/rotabase.txt
(rotabase only needed for FoldX 4.x; FoldX 5.x ignores it.)
"""
import argparse, subprocess, shutil, csv
from pathlib import Path
from statistics import mean, stdev

VARIANTS = ["Y1853C","A1843P","V1736A","M1783T","V1808A","V1665M",
            "R1751Q","L1664P","M1663K","P1806A","R1699L","R1699Q"]
CHAIN = "X"   # 1JNX chain; native BRCA1 numbering, offset 0
STRUCT = "1JNX_processed.pdb"
N_RUNS = 5

def run(cmd, cwd):
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd), timeout=1200)
    if r.returncode != 0:
        print("  FAILED rc=%d\n  STDERR: %s\n  STDOUT tail: %s" % (r.returncode, r.stderr[:400], r.stdout[-400:]))
    return r.returncode == 0

def parse_dif(dif):
    vals=[]
    for line in open(dif):
        line=line.strip()
        if not line or line.startswith(("Pdb","#")): continue
        p=line.split("\t")
        if len(p)>=2:
            try: vals.append(float(p[1]))
            except ValueError: pass
    return vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--foldx", required=True)
    ap.add_argument("--rotabase", default=None)
    ap.add_argument("--workdir", default="work")
    a=ap.parse_args()
    foldx=Path(a.foldx).resolve(); here=Path(__file__).parent.resolve()
    wd=(here/a.workdir); wd.mkdir(exist_ok=True)
    shutil.copy2(here/STRUCT, wd/STRUCT)
    if a.rotabase and Path(a.rotabase).exists():
        shutil.copy2(a.rotabase, wd/"rotabase.txt")

    # 1) RepairPDB once
    print("[1/2] RepairPDB on %s ..." % STRUCT)
    if not run([str(foldx),"--command=RepairPDB",f"--pdb={STRUCT}",f"--output-dir={wd}"], wd):
        print("RepairPDB failed; aborting."); return
    base=STRUCT.replace(".pdb","")
    repaired = wd/f"{base}_Repair.pdb"
    if not repaired.exists():
        cand=list(wd.glob("*_Repair.pdb"))
        repaired = cand[0] if cand else None
    if not repaired: print("No repaired PDB; aborting."); return
    rname=repaired.name; rbase=rname.replace(".pdb","")

    # 2) BuildModel per variant (n_runs=5) on the repaired structure
    rows=[]
    for v in VARIANTS:
        ref,pos,alt=v[0],v[1:-1],v[-1]
        mdir=wd/("bm_"+v); mdir.mkdir(exist_ok=True)
        shutil.copy2(repaired, mdir/rname)
        if a.rotabase and (wd/"rotabase.txt").exists():
            shutil.copy2(wd/"rotabase.txt", mdir/"rotabase.txt")
        (mdir/"individual_list.txt").write_text(f"{ref}{CHAIN}{pos}{alt};\n")
        print("[2/2] BuildModel %s ..." % v)
        ok=run([str(foldx),"--command=BuildModel",f"--pdb={rname}",
                "--mutant-file=individual_list.txt",f"--numberOfRuns={N_RUNS}",
                f"--output-dir={mdir}"], mdir)
        dif=mdir/f"Dif_{rbase}.fxout"
        if not dif.exists():
            c=list(mdir.glob("Dif_*.fxout")); dif=c[0] if c else None
        if ok and dif and dif.exists():
            vals=parse_dif(dif)
            m=round(mean(vals),4); sd=round(stdev(vals),4) if len(vals)>1 else 0.0
            rows.append((v,ref,pos,alt,m,sd,len(vals),";".join(map(str,vals))))
            print("      ddg_monomer = %.3f +/- %.3f (n=%d)" % (m,sd,len(vals)))
        else:
            rows.append((v,ref,pos,alt,"","",0,"")); print("      NO OUTPUT")

    out=here/"brct_foldx_ddg.csv"
    with open(out,"w",newline="") as f:
        w=csv.writer(f); w.writerow(["variant","ref_aa","position","alt_aa","foldx_ddg_monomer_mean","foldx_ddg_sd","n_runs","per_run_ddg"])
        w.writerows(rows)
    print("\nWrote %s (%d variants)" % (out,len(rows)))

if __name__=="__main__": main()
