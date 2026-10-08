#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 8 --mem 128gb --time=8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step3_n88.%j.log
# n88 copy of step3_cluster_composition.sh: cohort composition from
# results/clusters_{95,90,50}_n88.tsv (written by step25b; no createtsv, no MMseqs2).
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
DB=$WPS/catalog/db/LsPS_AA
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3

for T in 95 90 50; do
  TSV=$WPS/results/clusters_${T}_n88.tsv
  [ -s $TSV ] || { echo "missing $TSV"; exit 1; }
  echo "tier $T tsv lines: $(wc -l < $TSV)"
done

$PY << 'PYEOF'
import collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
wf22 = set(open(W+"/metadata/wf22_stems_44.txt").read().split())
wf23 = set(open(W+"/metadata/wf23_stems.txt").read().split())
def cohort(h):
    s = h.split("__")[0]
    return "WF22" if s in wf22 else ("WF23" if s in wf23 else "WF24")
for T in (95, 90, 50):
    comp = collections.defaultdict(set)
    size = collections.Counter()
    for line in open(f"{W}/results/clusters_{T}_n88.tsv"):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        comp[rep].add(cohort(mem)); size[rep] += 1
    pat = collections.Counter("+".join(sorted(v)) for v in comp.values())
    print(f"\n=== tier {T}: {len(comp)} clusters ===")
    for k, v in sorted(pat.items(), key=lambda x: -x[1]):
        print(f"  {k:20} {v:>10}  ({100.0*v/len(comp):5.2f}%)")
    singles = sum(1 for r in size if size[r] == 1)
    print(f"  singletons: {singles} ({100.0*singles/len(size):.2f}%)")
PYEOF
