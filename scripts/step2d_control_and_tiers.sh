#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 128gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step2d.%j.log
module load mmseqs2/13-45111
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
TMP=$SCRATCH/step2d_$SLURM_JOB_ID
mkdir -p $TMP
CPU=${SLURM_CPUS_ON_NODE:-32}

# ---- positive control: WF22 proteins (which BUILT the catalog) ----
CTRL=$WPS/catalog/wf22_control_200k.faa
if [ ! -s $CTRL ]; then
  for S in $(cat $WPS/metadata/wf22_stems.txt); do cat $V1/input_pep/$S.fasta; done | $PY -c "
import sys, random
random.seed(7)
p = 200000.0/18300000
keep=False; n=0
for line in sys.stdin:
    if line.startswith('>'):
        keep = random.random() < p
        if keep: n+=1
    if keep: sys.stdout.write(line)
sys.stderr.write('sampled %d\n' % n)
" > $CTRL
fi
echo "CONTROL query seqs: $(grep -c '^>' $CTRL)"

mmseqs easy-search $CTRL $V1/db/LsFMGC_AA_95_rep.fasta \
  $WPS/results/step2d_control_hits.m8 $TMP/c \
  --min-seq-id 0.95 -c 0.8 --cov-mode 1 --max-accept 1 -e 1e-5 -s 4 --threads $CPU

# ---- family-level: new proteins vs the 50% tier ----
mmseqs easy-search $WPS/catalog/subsample_200k.faa $V1/db/LsFMGC_AA_50_rep.fasta \
  $WPS/results/step2d_tier50_hits.m8 $TMP/t \
  --min-seq-id 0.5 -c 0.8 --cov-mode 1 --max-accept 1 -e 1e-5 -s 4 --threads $CPU

echo
echo "==================== RESULTS ===================="
CQ=$(grep -c '^>' $CTRL)
CH=$(cut -f1 $WPS/results/step2d_control_hits.m8 | sort -u | wc -l)
$PY -c "print('CONTROL (WF22 vs 95%% catalog): %d/%d = %.2f%%  [expect ~100%% if diagnostic is sound]' % ($CH,$CQ,100.0*$CH/$CQ))"

NQ=$(grep -c '^>' $WPS/catalog/subsample_200k.faa)
NH=$(cut -f1 $WPS/results/step2d_tier50_hits.m8 | sort -u | wc -l)
$PY -c "print('NEW vs 50%% tier (family level): %d/%d = %.2f%%' % ($NH,$NQ,100.0*$NH/$NQ))"
echo "NEW vs 95% tier (gene level): 58461/198842 = 29.40%  [from step2c]"

echo
echo "== length of mapped vs unmapped new proteins (95% search) =="
$PY << PYEOF
hits=set(l.split('\t')[0] for l in open("$WPS/results/step2c_hits.m8"))
import statistics
m=[];u=[]
name=None;L=0
for line in open("$WPS/catalog/subsample_200k.faa"):
    if line.startswith('>'):
        if name is not None:
            (m if name in hits else u).append(L)
        name=line[1:].strip(); L=0
    else:
        L+=len(line.strip())
if name is not None: (m if name in hits else u).append(L)
print("mapped:   n=%d median_aa=%d mean_aa=%d" % (len(m), statistics.median(m), sum(m)/len(m)))
print("unmapped: n=%d median_aa=%d mean_aa=%d" % (len(u), statistics.median(u), sum(u)/len(u)))
PYEOF
rm -rf $TMP
