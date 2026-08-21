#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 8 --mem 128gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step4.%j.log
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
wf22 = set(open(W+"/metadata/wf22_stems_44.txt").read().split())
wf23 = set(open(W+"/metadata/wf23_stems.txt").read().split())
def cohort(h):
    s = h.split("__")[0]
    return "WF22" if s in wf22 else ("WF23" if s in wf23 else "WF24")

for T in (95, 50):
    comp = collections.defaultdict(set)   # rep -> cohorts
    nsamp = collections.defaultdict(set)  # rep -> stems
    size = collections.Counter()
    for line in open(f"{W}/results/clusters_{T}.tsv"):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        comp[rep].add(cohort(mem))
        nsamp[rep].add(mem.split("__")[0])
        size[rep] += 1

    print(f"\n########## TIER {T} ##########")
    for label, keep in (("ALL", lambda r: True),
                        (">=2 members", lambda r: size[r] >= 2),
                        (">=2 samples", lambda r: len(nsamp[r]) >= 2),
                        (">=3 members AND >=2 samples",
                         lambda r: size[r] >= 3 and len(nsamp[r]) >= 2)):
        reps = [r for r in comp if keep(r)]
        pat = collections.Counter("+".join(sorted(comp[r])) for r in reps)
        n = len(reps)
        print(f"\n-- {label}: {n} clusters ({100.0*n/len(comp):.1f}% of all) --")
        for k, v in sorted(pat.items(), key=lambda x: -x[1]):
            print(f"   {k:18} {v:>9}  ({100.0*v/n:5.2f}%)")

    print("\n-- cluster size distribution --")
    sizes = sorted(size.values(), reverse=True)
    tot = sum(sizes)
    for frac in (0.001, 0.01, 0.1, 0.5):
        k = max(1, int(len(sizes)*frac))
        print(f"   top {frac*100:5.1f}% of clusters hold {100.0*sum(sizes[:k])/tot:5.2f}% of proteins")
    hist = collections.Counter(min(s, 11) for s in sizes)
    print("   size histogram (11 = 11+):")
    for s in sorted(hist):
        print(f"     {s if s<11 else '11+':>4}: {hist[s]}")
PYEOF
