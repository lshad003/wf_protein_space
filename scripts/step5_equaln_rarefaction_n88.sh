#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 8 --mem 200gb --time=8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step5_n88.%j.log
# n88 copy of step5_equaln_rarefaction.sh: results/clusters_95_n88.tsv, WF24
# stems without UHM586.41010 (35); same seeds, last pooled point 88.
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import collections, random
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
wf22 = [s for s in open(W+"/metadata/wf22_stems_44.txt").read().split()]
wf23 = [s for s in open(W+"/metadata/wf23_stems.txt").read().split()]
wf24 = [s for s in open(W+"/metadata/wf24_stems.txt").read().split() if s != "UHM586.41010"]
assert len(wf24) == 35
coh = {}
for s in wf22: coh[s]="WF22"
for s in wf23: coh[s]="WF23"
for s in wf24: coh[s]="WF24"

T = 95
print("reading cluster table, tier", T)
rep_id = {}
pairs = []            # (cluster_int, sample_str)
for line in open(f"{W}/results/clusters_{T}_n88.tsv"):
    rep, mem = line.rstrip("\n").split("\t")[:2]
    c = rep_id.get(rep)
    if c is None:
        c = len(rep_id); rep_id[rep] = c
    pairs.append((c, mem.split("__")[0]))
print("clusters:", len(rep_id), "member records:", len(pairs))

by_sample = collections.defaultdict(list)
for c, s in pairs:
    by_sample[s].append(c)

print("\n=== EQUAL-N COHORT COMPARISON (9 samples per cohort, 10 draws) ===")
res = collections.defaultdict(list)
for it in range(10):
    random.seed(100+it)
    sel = random.sample(wf22,9) + random.sample(wf23,9) + random.sample(wf24,9)
    selset = set(sel)
    cc = collections.defaultdict(set); cn = collections.defaultdict(set); cs = collections.Counter()
    for c, s in pairs:
        if s in selset:
            cc[c].add(coh[s]); cn[c].add(s); cs[c]+=1
    keep = [c for c in cc if cs[c]>=3 and len(cn[c])>=2]
    pat = collections.Counter("+".join(sorted(cc[c])) for c in keep)
    n = len(keep)
    for k,v in pat.items(): res[k].append(100.0*v/n)
    res["_n"].append(n)
print(f"clusters passing filter, mean: {sum(res['_n'])/10:.0f}")
for k in sorted(res, key=lambda x: -sum(res[x])):
    if k=="_n": continue
    v = res[k]
    print(f"  {k:18} mean {sum(v)/len(v):5.2f}%  (min {min(v):.2f}, max {max(v):.2f})")

print("\n=== RAREFACTION: distinct clusters vs samples added ===")
allsamples = wf22+wf23+wf24
for it in range(3):
    random.seed(200+it)
    order = random.sample(allsamples, len(allsamples))
    seen = set(); pts=[]
    for i, s in enumerate(order, 1):
        seen.update(by_sample[s])
        if i in (1,5,10,20,30,40,50,60,70,80,88): pts.append((i,len(seen)))
    print(" perm", it+1, " ".join(f"{i}:{n}" for i,n in pts))

print("\n=== RAREFACTION within each cohort (own samples only) ===")
for name, lst in (("WF22",wf22),("WF23",wf23),("WF24",wf24)):
    random.seed(300)
    order = random.sample(lst, len(lst))
    seen=set(); pts=[]
    for i, s in enumerate(order,1):
        seen.update(by_sample[s])
        if i in (1,3,5,9,20,30,36,44): pts.append((i,len(seen)))
    print(f" {name}: " + " ".join(f"{i}:{n}" for i,n in pts))
PYEOF
