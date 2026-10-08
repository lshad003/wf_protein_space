#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 8 --mem 200gb --time=8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step7b_n88.%j.log
# n88 copy of step7b_wf24_byrun.sh: reads results/clusters_95_n88.tsv; group sizes
# count only stems present in it, so UHM586.41010 is not counted.
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 << 'PYEOF'
import collections, numpy as np
W="/bigdata/stajichlab/lshad003/wf_protein_space"
run={}
for line in open("/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/samples.csv"):
    p=line.strip().split(",")
    if len(p)>=2 and p[0]!="METAGENOMEID": run[p[0]]=p[1].split("/")[0]
wf22=open(W+"/metadata/wf22_stems_44.txt").read().split()
wf23=open(W+"/metadata/wf23_stems.txt").read().split()
wf24=open(W+"/metadata/wf24_stems.txt").read().split()

rid={}; per=collections.defaultdict(list)
for line in open(W+"/results/clusters_95_n88.tsv"):
    rep,mem=line.rstrip("\n").split("\t")[:2]
    c=rid.get(rep)
    if c is None: c=len(rid); rid[rep]=c
    per[mem.split("__")[0]].append(c)
NC=len(rid); del rid
per={k:np.array(v,dtype=np.int32) for k,v in per.items()}

def rich(stems, target, seed=5):
    arrs=[per[s] for s in stems if s in per]
    tot=sum(len(a) for a in arrs)
    if tot<target: return None, tot
    rng=np.random.default_rng(seed)
    pool=np.concatenate(arrs); rng.shuffle(pool)
    seen=np.zeros(NC,dtype=bool); seen[pool[:target]]=True
    return int(seen.sum()), tot

groups={
 "WF22 (all 44)": wf22,
 "WF23 (all 9, run 0426)": wf23,
 "WF24 run 0426": [s for s in wf24 if run.get(s)=="UCB_20250426_M005990"],
 "WF24 run 0730": [s for s in wf24 if run.get(s)=="UCB_20250730_M006342"],
}
groups={k:[s for s in v if s in per] for k,v in groups.items()}
for T in (3_600_000, 8_000_000):
    print(f"\n=== richness at {T:,} proteins ===")
    for name,st in groups.items():
        r,tot=rich(st,T)
        print(f"  {name:26} n={len(st):>2} pool={tot:>12,}  " + (f"{r:>10,} clusters" if r else "insufficient"))
PYEOF
