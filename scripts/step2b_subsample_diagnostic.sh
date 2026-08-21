#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 32 --mem 96gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step2b_diag.%j.log
module load mmseqs2/13-45111
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
TMP=$SCRATCH/step2b_$SLURM_JOB_ID
mkdir -p $TMP $WPS/results
CPU=${SLURM_CPUS_ON_NODE:-32}

QRY=$WPS/catalog/subsample_200k.faa
if [ ! -s $QRY ]; then
  cat $WPS/catalog/input_pep/*.fasta | $PY -c "
import sys, random
random.seed(42)
p = 200000.0/24151134
keep=False; n=0
for line in sys.stdin:
    if line.startswith('>'):
        keep = random.random() < p
        if keep: n+=1
    if keep: sys.stdout.write(line)
sys.stderr.write('sampled %d\n' % n)
" > $QRY
fi
echo "query seqs: $(grep -c '^>' $QRY)"

mmseqs easy-search $QRY $V1/db/LsFMGC_AA_95_rep.fasta \
  $WPS/results/step2b_hits.m8 $TMP \
  --min-seq-id 0.95 -c 0.8 --cov-mode 1 \
  --max-accept 1 -e 1e-5 --threads $CPU

echo "== RESULT =="
Q=$(grep -c '^>' $QRY)
H=$(cut -f1 $WPS/results/step2b_hits.m8 | sort -u | wc -l)
echo "queries: $Q ; with hit: $H"
$PY -c "print('mapping rate: %.2f%%' % (100.0*$H/$Q))"

echo "== per-year =="
grep '^>' $QRY | sed 's/>//; s/__.*//' | sort | uniq -c > $TMP/q_by_stem.txt
cut -f1 $WPS/results/step2b_hits.m8 | sed 's/__.*//' | sort -u -o /dev/null 2>/dev/null
cut -f1 $WPS/results/step2b_hits.m8 | sort -u | sed 's/__.*//' | sort | uniq -c > $TMP/h_by_stem.txt
$PY - << 'PYEOF'
import os
tmp=os.environ.get('TMP') or os.path.expandvars('$SCRATCH')
def load(f):
    d={}
    for line in open(f):
        c,s=line.split(); d[s]=int(c)
    return d
PYEOF
paste <(sort -k2 $TMP/q_by_stem.txt) <(sort -k2 $TMP/h_by_stem.txt) | head -50
rm -rf $TMP
