#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 4 --mem 16gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step6b.%j.log
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
FEC="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal"
wf22=set(open(W+"/metadata/wf22_stems_44.txt").read().split())
wf23=set(open(W+"/metadata/wf23_stems.txt").read().split())
wf24=set(open(W+"/metadata/wf24_stems.txt").read().split())

run={}
for line in open(FEC+"/samples.csv"):
    p=line.strip().split(",")
    if len(p)<2 or p[0]=="METAGENOMEID": continue
    run[p[0]]=p[1].split("/")[0]

prot={}
for line in open(W+"/results/step1_protein_counts_45.tsv"):
    f=line.split("\t")
    if f[0]!="stem": prot[f[0]]=int(f[3])
V1PEP="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1/input_pep"
for s in wf22:
    fp=f"{V1PEP}/{s}.fasta"
    if os.path.exists(fp):
        prot[s]=sum(1 for l in open(fp) if l.startswith(">"))

def lhist(path):
    n=0; wsum=0
    if not os.path.exists(path): return None,None
    for line in open(path):
        if line.startswith("#"): continue
        a=line.split()
        if len(a)>=2:
            L=int(a[0]); c=int(a[1]); n+=c; wsum+=L*c
    return n, (wsum/n if n else 0)

def gc(path):
    if not os.path.exists(path): return None
    for line in open(path):
        if line.startswith("#Mean"): return float(line.split()[1])
    return None

rows=[]
for s in sorted(set(list(wf22)+list(wf23)+list(wf24))):
    coh = "WF22" if s in wf22 else "WF23" if s in wf23 else "WF24"
    base=f"{FEC}/results/{s}"
    raw_n, raw_len = lhist(f"{base}/raw_reads_stats/{s}_R/lhist.txt")
    cln_n, cln_len = lhist(f"{base}/clean_reads_stats/{s}_R/lhist.txt")
    g = gc(f"{base}/clean_reads_stats/{s}_R/gchist.txt")
    r = run.get(s,"NA")
    center = r.split("_")[0] if r!="NA" else "NA"
    plat = "AVITI" if "AVITI" in r else "Illumina"
    rows.append(dict(stem=s, cohort=coh, run=r, center=center, platform=plat,
        raw_reads=raw_n, clean_reads=cln_n, pct_kept=(100.0*cln_n/raw_n if raw_n else None),
        clean_len_mean=cln_len, gc_mean=g, n_proteins=prot.get(s)))

out=W+"/results/batch_table_89.tsv"
cols=["stem","cohort","run","center","platform","raw_reads","clean_reads","pct_kept","clean_len_mean","gc_mean","n_proteins"]
with open(out,"w") as fh:
    fh.write("\t".join(cols)+"\n")
    for r in rows:
        fh.write("\t".join("NA" if r[c] is None else (f"{r[c]:.2f}" if isinstance(r[c],float) else str(r[c])) for c in cols)+"\n")
print("wrote",out,"rows",len(rows))
print("missing clean_reads:", sum(1 for r in rows if r["clean_reads"] is None))

def summarize(key):
    g=collections.defaultdict(list)
    for r in rows:
        if r["clean_reads"]: g[r[key]].append(r)
    print(f"\n=== by {key} ===")
    print(f"{key:28} {'n':>3} {'clean_reads_median':>19} {'len':>6} {'GC':>6} {'proteins_median':>16}")
    for k in sorted(g):
        v=g[k]; med=lambda f: sorted(x[f] for x in v if x[f])[len(v)//2]
        print(f"{str(k):28} {len(v):>3} {med('clean_reads'):>19,} {med('clean_len_mean'):>6.1f} {med('gc_mean'):>6.2f} {med('n_proteins'):>16,}")
summarize("cohort"); summarize("run"); summarize("platform")

print("\n=== KEY CONTRAST 1: WF23 vs WF24 within UCB_20250426_M005990 ===")
sub=[r for r in rows if r["run"]=="UCB_20250426_M005990" and r["clean_reads"]]
for c in ("WF23","WF24"):
    v=[r for r in sub if r["cohort"]==c]
    if v: print(f"  {c}: n={len(v)} median_clean={sorted(x['clean_reads'] for x in v)[len(v)//2]:,} median_proteins={sorted(x['n_proteins'] for x in v)[len(v)//2]:,}")
print("\n=== KEY CONTRAST 2: WF24 across its two runs (batch effect estimate) ===")
for rn in ("UCB_20250426_M005990","UCB_20250730_M006342"):
    v=[r for r in rows if r["cohort"]=="WF24" and r["run"]==rn and r["clean_reads"]]
    if v: print(f"  {rn}: n={len(v)} median_clean={sorted(x['clean_reads'] for x in v)[len(v)//2]:,} median_proteins={sorted(x['n_proteins'] for x in v)[len(v)//2]:,}")
PYEOF
