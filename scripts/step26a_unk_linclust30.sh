#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 256gb --time=24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step26a.%j.log
# Step 26a: metagRoot/Aplakidou 2026 comparable families, a comparability
# analysis only. The 50% families (results/unk_clusters_50_n88.tsv) stay frozen.
# Input: results/supported_U_reps_n88.faa (585,347 unknown representatives).
# Filters reproducible here:
#   F1 complete ORF, or partial but not within 10 nt of a contig end
#      (Prodigal partial= flags, read from the per-sample aa.fa.gz as in Step 9;
#      for left-edge partials the Prodigal start coordinate is also checked)
#   F2 length >= 35 aa (terminal * not counted)
# Clustering: mmseqs linclust --min-seq-id 0.3 -c 0.8 --cov-mode 0.
# The 40% run (step14d) crashed with "Sequence db size != result db size" using
# a tmp dir on shared storage; here tmp is node-local $SCRATCH (workspace/scratch)
# and memory is 256 GB. First attempt, job 29658478, failed: /scratch/$USER not writable.
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
IN=$W/results/supported_U_reps_n88.faa
FAA=$W/results/unk_filtered_aplakidou.faa
TSV=$W/results/unk_clusters_30.tsv
SIZES=$W/results/unk_family_sizes_30.tsv
for f in $FAA $TSV $SIZES; do
  if [ -e $f ]; then echo "REFUSING: $f exists"; exit 1; fi
done

$PY << 'PYEOF'
import gzip, os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
IN=W+"/results/supported_U_reps_n88.faa"
OUT=W+"/results/unk_filtered_aplakidou.faa"

seqs={}; order=[]; cur=None
for line in open(IN):
    if line[0]==">":
        cur=line[1:].split()[0]; order.append(cur); seqs[cur]=[]
    else:
        seqs[cur].append(line.strip())
seqs={k:"".join(v) for k,v in seqs.items()}
print("input reps:", len(order), flush=True)
stems=sorted({k.split("__",1)[0] for k in order})
print("stems with reps:", len(stems), flush=True)

info={}   # id -> (partial, start, end)
for i,s in enumerate(stems,1):
    p = f"{NEW}/{s}.aa.fa.gz" if os.path.exists(f"{NEW}/{s}.aa.fa.gz") else f"{OLD}/{s}.aa.fa.gz"
    with gzip.open(p,"rt") as fh:
        for line in fh:
            if line[0]!=">": continue
            f=line[1:].split(" # ")
            pid=s+"__"+f[0].split()[0]
            if pid not in seqs: continue
            j=line.find("partial=")
            code=line[j+8:j+10] if j>=0 else "NA"
            info[pid]=(code,int(f[1]),int(f[2]))
    if i%10==0: print(" ",i,"stems done", flush=True)
missing=[k for k in order if k not in info]
print("reps without a Prodigal header:", len(missing), flush=True)
if missing: raise SystemExit("ABORT: missing flags")

pc=collections.Counter(v[0] for v in info.values())
print("\npartial flags among input reps:")
for k in sorted(pc):
    lab={"00":"complete","10":"no start (left edge)","01":"no stop (right edge)","11":"both edges"}.get(k,k)
    print(f"  {k} {lab:22} {pc[k]:>9,}  {100.0*pc[k]/len(order):6.2f}%")
# Prodigal only flags a gene partial when it runs off the contig, so a partial
# gene sits at the edge. Check the left edge from coordinates (start <= 10).
left=[v for v in info.values() if v[0][0]=="1"]
far=sum(1 for v in left if v[1]>10)
print(f"left-edge partials with start > 10 nt: {far:,} of {len(left):,}")
print("right-edge distance needs contig lengths; not checked (Prodigal partial= means the ORF runs off that end)")

def L(k): return len(seqs[k].rstrip("*"))
keep1=lambda k: info[k][0]=="00" or (info[k][0][0]=="1" and info[k][1]>10 and info[k][0][1]=="0")
keep2=lambda k: L(k)>=35
n=len(order)
f1=sum(1 for k in order if not keep1(k))
f2=sum(1 for k in order if not keep2(k))
both=sum(1 for k in order if not keep1(k) and not keep2(k))
f2_after_f1=sum(1 for k in order if keep1(k) and not keep2(k))
kept=[k for k in order if keep1(k) and keep2(k)]
print(f"\nF1 complete / not within 10 nt of edge: removes {f1:,} ({100.0*f1/n:.2f}%)")
print(f"F2 length >= 35 aa:                   removes {f2:,} ({100.0*f2/n:.2f}%) on its own")
print(f"   removed by both filters: {both:,}; F2 removes {f2_after_f1:,} after F1")
print(f"kept: {len(kept):,} of {n:,} ({100.0*len(kept)/n:.2f}%)")
with open(OUT,"w") as fh:
    for k in kept: fh.write(f">{k}\n{seqs[k]}\n")
print("wrote", OUT, flush=True)
PYEOF
if [ ! -s $FAA ]; then echo "ABORT: no filtered fasta"; exit 1; fi

module load mmseqs2/13-45111
mmseqs version
module load workspace/scratch
T=$SCRATCH/step26a
mkdir -p $T
cp $FAA $T/in.faa
mmseqs createdb $T/in.faa $T/U30
mmseqs linclust $T/U30 $T/U30_clu $T/tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 \
  --threads $SLURM_CPUS_PER_TASK --split-memory-limit 200G
mmseqs createtsv $T/U30 $T/U30 $T/U30_clu $T/clusters_30.tsv --threads $SLURM_CPUS_PER_TASK
if [ ! -s $T/clusters_30.tsv ]; then echo "ABORT: linclust produced no tsv"; exit 1; fi
cp $T/clusters_30.tsv $TSV

echo "members in tsv: $(wc -l < $TSV)  filtered seqs: $(grep -c '>' $FAA)"
awk -F'\t' '{n[$1]++} END{
  for(r in n){t++; m+=n[r]; s=n[r]; if(s>=3){a++;pa+=s} if(s>=25){b++;pb+=s} if(s>=100){d++;pd+=s}}
  printf "families total\t%d\n", t
  printf ">=3 members\t%d\tproteins %d (%.2f%%)\n", a,pa,100*pa/m
  printf ">=25 members\t%d\tproteins %d (%.2f%%)\n", b,pb,100*pb/m
  printf ">=100 members\t%d\tproteins %d (%.2f%%)\n", d,pd,100*pd/m
}' $TSV > $SIZES
echo "linclust 30% id, -c 0.8 --cov-mode 0:"
cat $SIZES
rm -rf $T
echo "DONE step26a"
