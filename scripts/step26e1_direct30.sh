#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 128gb --time=24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step26e1.%j.log
# Step 26e1: direct 30% clustering of member proteins (gold standard for family size).
# Input: every member of the 384,260 Aplakidou-filtered unknown reps' 95% clusters
# (results/clusters_95_n88.tsv) that is complete (partial=00) and >= 35 aa, sequences
# from the per-sample Prodigal aa.fa.gz (same source and filter as step26d2).
# Clustering: mmseqs linclust --min-seq-id 0.3 -c 0.8 --cov-mode 0 (as step26a and the
# Aplakidou paper; --cov-mode 0 = 80% bidirectional, deliberately not the search default).
# Outputs results/unk_direct30_{clusters,sizes,summary}.tsv; nothing existing is overwritten.
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
for f in $W/results/unk_direct30_clusters.tsv $W/results/unk_direct30_sizes.tsv $W/results/unk_direct30_summary.tsv; do
  if [ -e $f ]; then echo "REFUSING: $f exists"; exit 1; fi
done
module load workspace/scratch
T=$SCRATCH/step26e1
mkdir -p $T

$PY - $T/members.faa << 'PYEOF'
import gzip, os, sys
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
reps={l[1:].split()[0] for l in open(W+"/results/unk_filtered_aplakidou.faa") if l[0]==">"}
mem=set()
for line in open(W+"/results/clusters_95_n88.tsv"):
    r,m=line.rstrip("\n").split("\t")
    if r in reps: mem.add(m)
print("filtered reps:", len(reps), "| member proteins:", len(mem), flush=True)
def fasta(fh):
    h=None; s=[]
    for line in fh:
        if line[0]==">":
            if h is not None: yield h,"".join(s)
            h=line[1:].rstrip("\n"); s=[]
        else: s.append(line.strip())
    if h is not None: yield h,"".join(s)
seen=0; kept=0
out=open(sys.argv[1],"w")
for s in sorted({m.split("__",1)[0] for m in mem}):
    p = f"{NEW}/{s}.aa.fa.gz" if os.path.exists(f"{NEW}/{s}.aa.fa.gz") else f"{OLD}/{s}.aa.fa.gz"
    with gzip.open(p,"rt") as fh:
        for h,seq in fasta(fh):
            pid=s+"__"+h.split(" # ")[0].split()[0]
            if pid not in mem: continue
            seen+=1
            j=h.find("partial="); code=h[j+8:j+10] if j>=0 else "NA"
            seq=seq.rstrip("*")
            if code=="00" and len(seq)>=35:
                kept+=1; out.write(f">{pid}\n{seq}\n")
out.close()
print("members found:", seen, "| complete >=35 aa written:", kept, flush=True)
if seen!=len(mem): raise SystemExit("ABORT: members missing from Prodigal files")
PYEOF
if [ ! -s $T/members.faa ]; then echo "ABORT: no member fasta"; exit 1; fi

module load mmseqs2/13-45111
mmseqs version
mmseqs createdb $T/members.faa $T/M
mmseqs linclust $T/M $T/M_clu $T/tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 \
  --threads $SLURM_CPUS_PER_TASK --split-memory-limit 100G
mmseqs createtsv $T/M $T/M $T/M_clu $T/clusters.tsv --threads $SLURM_CPUS_PER_TASK
if [ ! -s $T/clusters.tsv ]; then echo "ABORT: linclust produced no tsv"; exit 1; fi
cp $T/clusters.tsv $W/results/unk_direct30_clusters.tsv

$PY - $T/members.faa << 'PYEOF'
import sys, hashlib, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
h={}; k=None
for l in open(sys.argv[1]):
    if l[0]==">": k=l[1:].strip()
    else: h[k]=hashlib.md5(l.strip().encode()).digest()
n=collections.Counter(); u=collections.defaultdict(set)
for l in open(W+"/results/unk_direct30_clusters.tsv"):
    f,m=l.rstrip("\n").split("\t"); n[f]+=1; u[f].add(h[m])
nu={f:len(v) for f,v in u.items()}
print("proteins clustered:", sum(n.values()), "of", len(h), "| families:", len(n))
with open(W+"/results/unk_direct30_sizes.tsv","w") as fh:
    fh.write("family_rep\tn_proteins\tn_unique\n")
    for f in sorted(n, key=lambda x:(-n[x],x)): fh.write(f"{f}\t{n[f]}\t{nu[f]}\n")
rows=[("count","families_total",">=3",">=25",">=100")]
for lab,c in (("protein records",n),("exact-unique sequences",nu)):
    rows.append((lab,len(n))+tuple(sum(1 for f in n if c[f]>=t) for t in (3,25,100)))
rows.append(("step26a2 proteins via 95% reps (complete >=35 aa)",358556,271708,14022,236))
rows.append(("step26d2 unique via 95% reps (complete >=35 aa)",358556,"NA",45,1))
with open(W+"/results/unk_direct30_summary.tsv","w") as fh:
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
for r in rows: print("\t".join(map(str,r)))
PYEOF
rm -rf $T
echo "DONE step26e1"
