#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 64 --mem 400gb --time=3-00:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step26e2.%j.log
# Step 26e2: positive control for the step26b3 NMPFams2 search.
# 5,000 random NMPFams2 metag members that are not family representatives
# (reservoir sample, seed 26) searched against the same representative database
# with exactly the step26b3 settings: mmseqs 13-45111 search -e 1e-3
# --min-seq-id 0.3 -c 0.7 --cov-mode 0, default -s 5.7.
# Recovered = hit to any representative; correct = hit to its own family's representative.
# Outputs results/nmpfams2_posctrl_{queries.faa,hits.tsv,summary.tsv}.
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
BZ=$W/catalog/nmpfams2/metag_clusters25_names.tsv.bz2
for f in $W/results/nmpfams2_posctrl_queries.faa $W/results/nmpfams2_posctrl_hits.tsv $W/results/nmpfams2_posctrl_summary.tsv; do
  if [ -e $f ]; then echo "REFUSING: $f exists"; exit 1; fi
done
module load workspace/scratch
T=$SCRATCH/step26e2
mkdir -p $T

bzcat $BZ | $PY -c "
import sys, random
random.seed(26)
reps=open('$T/reps.faa','w'); res=[]; n=0; nrep=0; seen=set()
for line in sys.stdin:
    f=line.rstrip('\n').split('\t')
    if len(f)<4: continue
    fam,rep,mem,seq=f[0],f[1],f[2],f[3]
    if rep==mem:
        if fam in seen: continue
        seen.add(fam); nrep+=1; reps.write(f'>{fam}\n{seq}\n'); continue
    n+=1
    if len(res)<5000: res.append((fam,mem,seq))
    else:
        j=random.randrange(n)
        if j<5000: res[j]=(fam,mem,seq)
reps.close()
with open('$T/q.faa','w') as fh:
    for i,(fam,mem,seq) in enumerate(res): fh.write(f'>q{i+1}|{fam}\n{seq}\n')
print(f'representatives\t{nrep}'); print(f'non-representative members\t{n}'); print(f'sampled\t{len(res)}')
"
if [ ! -s $T/q.faa ]; then echo "ABORT: no control queries"; exit 1; fi
cp $T/q.faa $W/results/nmpfams2_posctrl_queries.faa

module load mmseqs2/13-45111
mmseqs version
mmseqs createdb $T/q.faa $T/Q
mmseqs createdb $T/reps.faa $T/R
mmseqs search $T/Q $T/R $T/QR $T/tmp -e 1e-3 --min-seq-id 0.3 -c 0.7 --cov-mode 0 \
  --threads $SLURM_CPUS_PER_TASK --split-memory-limit 300G
mmseqs convertalis $T/Q $T/R $T/QR $T/QR.tsv \
  --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits --threads $SLURM_CPUS_PER_TASK
if [ ! -e $T/QR.tsv ]; then echo "ABORT: no search output"; exit 1; fi
cp $T/QR.tsv $W/results/nmpfams2_posctrl_hits.tsv

$PY -c "
import collections
W='$W'
q=[l[1:].split()[0] for l in open(W+'/results/nmpfams2_posctrl_queries.faa') if l[0]=='>']
L={}; k=None
for l in open(W+'/results/nmpfams2_posctrl_queries.faa'):
    if l[0]=='>': k=l[1:].split()[0]
    else: L[k]=len(l.strip().rstrip('*'))
hits=collections.defaultdict(list)
for l in open(W+'/results/nmpfams2_posctrl_hits.tsv'):
    f=l.split('\t'); hits[f[0]].append(f[1])
anyhit=sum(1 for x in q if hits[x]); own=sum(1 for x in q if x.split('|',1)[1] in hits[x])
top=sum(1 for x in q if hits[x] and hits[x][0]==x.split('|',1)[1])
q35=[x for x in q if L[x]>=35]
own35=sum(1 for x in q35 if x.split('|',1)[1] in hits[x])
rows=[('measure','n','of','pct'),
 ('hit any representative',anyhit,len(q),f'{100.0*anyhit/len(q):.2f}'),
 ('hit own family representative',own,len(q),f'{100.0*own/len(q):.2f}'),
 ('top hit is own family',top,len(q),f'{100.0*top/len(q):.2f}'),
 ('own family, queries >=35 aa',own35,len(q35),f'{100.0*own35/len(q35):.2f}')]
with open(W+'/results/nmpfams2_posctrl_summary.tsv','w') as fh:
    for r in rows: fh.write('\t'.join(map(str,r))+'\n')
for r in rows: print('\t'.join(map(str,r)))
"
rm -rf $T
echo "DONE step26e2"
