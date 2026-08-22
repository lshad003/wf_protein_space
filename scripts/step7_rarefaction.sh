#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 8 --mem 200gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step7.%j.log
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import collections, random
import numpy as np
W="/bigdata/stajichlab/lshad003/wf_protein_space"
wf22=open(W+"/metadata/wf22_stems_44.txt").read().split()
wf23=open(W+"/metadata/wf23_stems.txt").read().split()
wf24=open(W+"/metadata/wf24_stems.txt").read().split()
coh={}
for s in wf22: coh[s]="WF22"
for s in wf23: coh[s]="WF23"
for s in wf24: coh[s]="WF24"

print("reading clusters_95.tsv ...", flush=True)
rid={}
per_sample=collections.defaultdict(list)
for line in open(W+"/results/clusters_95.tsv"):
    rep,mem=line.rstrip("\n").split("\t")[:2]
    c=rid.get(rep)
    if c is None:
        c=len(rid); rid[rep]=c
    per_sample[mem.split("__")[0]].append(c)
NC=len(rid)
print("clusters:",NC,"samples:",len(per_sample), flush=True)
per_sample={k:np.array(v,dtype=np.int32) for k,v in per_sample.items()}
del rid

def curve(arrays, seed, points):
    """pool arrays, shuffle, report distinct clusters at each protein count"""
    rng=np.random.default_rng(seed)
    pool=np.concatenate(arrays)
    rng.shuffle(pool)
    seen=np.zeros(NC,dtype=bool)
    out=[]; n=0; idx=0
    for p in points:
        if p>len(pool): break
        chunk=pool[idx:p]
        seen[chunk]=True
        idx=p
        out.append((p,int(seen.sum())))
    return out

POINTS=[100_000,250_000,500_000,1_000_000,2_000_000,3_000_000,3_600_000,
        5_000_000,7_500_000,10_000_000,15_000_000,16_000_000,20_000_000]

print("\n=== RAREFACTION BY PROTEINS SAMPLED (per cohort) ===", flush=True)
rows=[]
for name,lst in (("WF22",wf22),("WF23",wf23),("WF24",wf24)):
    arrs=[per_sample[s] for s in lst if s in per_sample]
    tot=sum(len(a) for a in arrs)
    print(f"\n{name}: {len(arrs)} samples, {tot:,} proteins")
    reps=[]
    for it in range(3):
        reps.append(dict(curve(arrs, 1000+it, POINTS)))
    for p in POINTS:
        vals=[r[p] for r in reps if p in r]
        if vals:
            m=sum(vals)/len(vals)
            print(f"   {p:>12,} proteins -> {m:>12,.0f} clusters")
            rows.append((name,p,m,min(vals),max(vals)))

print("\n=== FAIR COMPARISON at 3,600,000 proteins (WF23's total) ===")
for name,p,m,lo,hi in rows:
    if p==3_600_000: print(f"   {name}: {m:,.0f} clusters  (range {lo:,}-{hi:,})")

print("\n=== POOLED, all 89 samples ===")
arrs=[per_sample[s] for s in per_sample]
tot=sum(len(a) for a in arrs)
print(f"total proteins: {tot:,}")
BIG=POINTS+[25_000_000,30_000_000,35_000_000,40_000_000,40_550_595]
pooled=dict(curve(arrs,7,BIG))
prev=None
for p in BIG:
    if p in pooled:
        d=""
        if prev: d=f"   (+{pooled[p]-prev[1]:,} over last {p-prev[0]:,})"
        print(f"   {p:>12,} proteins -> {pooled[p]:>12,} clusters{d}")
        prev=(p,pooled[p])
        rows.append(("POOLED",p,pooled[p],pooled[p],pooled[p]))

out=W+"/results/rarefaction_95.tsv"
with open(out,"w") as fh:
    fh.write("group\tproteins_sampled\tclusters_mean\tclusters_min\tclusters_max\n")
    for r in rows: fh.write("\t".join(str(x) for x in r)+"\n")
print("\nwrote",out)
PYEOF
