#!/usr/bin/bash -l
#SBATCH -J step34c_strict_eu
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 128gb --time=24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34c_strict_eu.%j.log
# Step 34c: 30% families of the 41,392 strict EU reps, recurrence, NMPFams2, Cse2 control.
# Clustering exactly as step26e1: mmseqs linclust --min-seq-id 0.3 -c 0.8 --cov-mode 0
#   (--cov-mode 0 = 80% bidirectional, as step26a/step26e1 and Aplakidou; deliberately not the search default),
#   --split-memory-limit 100G, other parameters at MMseqs2 13-45111 defaults.
# Difference from step26e1: the units clustered are the strict EU 95% reps themselves (sequence from
#   catalog/db/LsPS_AA_95_rep.fasta, terminal "*" stripped), not their complete >= 35 aa member proteins.
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
for f in $W/results/step34c_strict_eu_families.tsv $W/results/step34c_strict_eu_nmpfams2.tsv $W/results/step34c_strict_eu_clusters30.tsv; do
  if [ -e $f ]; then echo "REFUSING: $f exists"; exit 1; fi
done
hostname; date
echo "step26e1 parameters (from $W/scripts/step26e1_direct30.sh):"
grep -E "^mmseqs linclust|^  --threads" $W/scripts/step26e1_direct30.sh
module load workspace/scratch
T=$SCRATCH/step34c
mkdir -p $T
$PY $W/scripts/step34c_prep.py $T/strict_eu.faa || { echo "exit 1 (prep)"; exit 1; }
module load mmseqs2/13-45111
mmseqs version
mmseqs createdb $T/strict_eu.faa $T/M
mmseqs linclust $T/M $T/M_clu $T/tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 \
  --threads $SLURM_CPUS_PER_TASK --split-memory-limit 100G
mmseqs createtsv $T/M $T/M $T/M_clu $T/clusters.tsv --threads $SLURM_CPUS_PER_TASK
if [ ! -s $T/clusters.tsv ]; then echo "ABORT: linclust produced no tsv"; exit 1; fi
cp $T/clusters.tsv $W/results/step34c_strict_eu_clusters30.tsv
$PY $W/scripts/step34c_families.py $W/results/step34c_strict_eu_clusters30.tsv
echo "exit $?"
rm -rf $T
date
