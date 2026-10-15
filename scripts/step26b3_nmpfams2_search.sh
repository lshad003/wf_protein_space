#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 64 --mem 400gb --time=3-00:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step26b3.%x.%j.log
# Step 26b3: search the Aplakidou-filtered unknowns against NMPFams2 families.
# Usage: sbatch -J nmp_metag scripts/step26b3_nmpfams2_search.sh metag
#        sbatch -J nmp_iso   scripts/step26b3_nmpfams2_search.sh iso
# Source file catalog/nmpfams2/<src>_clusters25_names.tsv.bz2, columns:
#   family, rep ID, member ID, sequence; representative rows have col2 == col3.
# Queries: results/unk_filtered_aplakidou.faa (384,260).
# Search: mmseqs search -e 1e-3 --min-seq-id 0.3 -c 0.7 --cov-mode 0 (default -s 5.7).
#   --cov-mode 0 (bidirectional 70%) deviates from the CLAUDE.md --cov-mode 2 default
#   on purpose: NMPFams2 families are built at 80% bidirectional coverage, so
#   family membership needs coverage of both query and target. Agreed 2026-10-09.
# Sensitivity: a random 10,000 queries (seed 26) searched against ALL members.
# No _n88 file and no README.md is touched.
SRC=$1
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
D=$W/catalog/nmpfams2
BZ=$D/${SRC}_clusters25_names.tsv.bz2
Q=$W/results/unk_filtered_aplakidou.faa
HITS=$W/results/nmpfams2_${SRC}_rep_hits.tsv
SHITS=$W/results/nmpfams2_${SRC}_all_hits_10k.tsv
SUB=$W/results/nmpfams2_10k_queries.txt
case $SRC in
  metag) WANT=677b976daebebacc289667d9d1fd48b6 ;;
  iso)   WANT=d6ae4851bdb76667b233bb4a3bb2dc19 ;;
  *) echo "ABORT: source must be metag or iso"; exit 1 ;;
esac
for f in $HITS $SHITS; do
  if [ -e $f ]; then echo "REFUSING: $f exists"; exit 1; fi
done
if [ "$SRC" = "iso" ]; then
  grep -aq "^$WANT " $W/logs/step26b2.log || { echo "ABORT: logs/step26b2.log has no md5 match for $WANT"; exit 1; }
fi
GOT=$(md5sum $BZ | cut -d' ' -f1)
echo "md5 $GOT want $WANT"
if [ "$GOT" != "$WANT" ]; then echo "ABORT: md5 mismatch"; exit 1; fi

module load mmseqs2/13-45111
mmseqs version
module load workspace/scratch
T=$SCRATCH/step26b3_$SRC
mkdir -p $T
cp $Q $T/q.faa

# 1. representatives (one per family) and all members, one pass
bzcat $BZ | $PY -c "
import sys
reps=open('$T/reps.faa','w'); allm=open('$T/all.faa','w')
fams=set(); nrep=0; nmem=0; repfam=set(); dup=0
for line in sys.stdin:
    f=line.rstrip('\n').split('\t')
    if len(f)<4: continue
    fam,rep,mem,seq=f[0],f[1],f[2],f[3]
    fams.add(fam); nmem+=1
    allm.write(f'>{fam}|m{nmem}\n{seq}\n')
    if rep==mem:
        if fam in repfam: dup+=1; continue
        repfam.add(fam); nrep+=1
        reps.write(f'>{fam}\n{seq}\n')
reps.close(); allm.close()
print(f'families\t{len(fams)}'); print(f'members\t{nmem}')
print(f'representatives\t{nrep}'); print(f'families without a rep row\t{len(fams)-len(repfam)}')
print(f'extra rep rows in a family (skipped)\t{dup}')
" | tee $T/extract_counts.tsv
if [ ! -s $T/reps.faa ]; then echo "ABORT: no representatives extracted"; exit 1; fi
cp $T/extract_counts.tsv $W/results/nmpfams2_${SRC}_extract_counts.tsv

# 2. all queries vs representatives
mmseqs createdb $T/q.faa $T/Q
mmseqs createdb $T/reps.faa $T/R
mmseqs search $T/Q $T/R $T/QR $T/tmp1 -e 1e-3 --min-seq-id 0.3 -c 0.7 --cov-mode 0 \
  --threads $SLURM_CPUS_PER_TASK --split-memory-limit 300G
mmseqs convertalis $T/Q $T/R $T/QR $T/QR.tsv \
  --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits --threads $SLURM_CPUS_PER_TASK
if [ ! -e $T/QR.tsv ]; then echo "ABORT: no rep search output"; exit 1; fi
cp $T/QR.tsv $HITS
echo "rep hits rows: $(wc -l < $HITS); queries with a hit: $(cut -f1 $HITS | sort -u | wc -l)"

# 3. sensitivity: random 10,000 queries vs all members (same subset for every source)
if [ ! -s $SUB ]; then
  $PY -c "
import random
ids=[l[1:].split()[0] for l in open('$Q') if l[0]=='>']
random.seed(26); s=random.sample(ids,10000)
open('$SUB','w').write('\n'.join(s)+'\n')"
fi
$PY -c "
keep=set(l.strip() for l in open('$SUB'))
out=open('$T/q10k.faa','w'); w=False
for l in open('$T/q.faa'):
    if l[0]=='>': w=l[1:].split()[0] in keep
    if w: out.write(l)"
echo "10k subset seqs: $(grep -c '>' $T/q10k.faa)"
mmseqs createdb $T/q10k.faa $T/Q10
mmseqs createdb $T/all.faa $T/A
mmseqs search $T/Q10 $T/A $T/QA $T/tmp2 -e 1e-3 --min-seq-id 0.3 -c 0.7 --cov-mode 0 \
  --threads $SLURM_CPUS_PER_TASK --split-memory-limit 300G
mmseqs convertalis $T/Q10 $T/A $T/QA $T/QA.tsv \
  --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits --threads $SLURM_CPUS_PER_TASK
if [ ! -e $T/QA.tsv ]; then echo "ABORT: no all-member search output"; exit 1; fi
cp $T/QA.tsv $SHITS
rm -rf $T
echo "DONE step26b3 $SRC"
