#!/bin/bash
#SBATCH -J step33f_classify
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33f_classify.%j.log
# Step 33f classify: collect 18S hits >= 500 bp (scripts/step33f_collect.py), top hit against SILVA SSU 138.1 NR99
# (/srv/projects/db/SILVA/SSU/138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta, 510,508 sequences, DNA) with
# vsearch 2.21.1 --usearch_global --id 0.70 --strand both --maxaccepts 20 --maxrejects 100 --top_hits_only
# --notrunclabels --output_no_hits (database staged to node-local $SCRATCH), then scripts/step33f_summary.py.
# Outputs (refuse to overwrite): results/step33f_18S_ge500.fa, results/step33f_18S_silva.tsv (query, target with
# SILVA taxonomy, id, alnlen, ql, qcov), results/step33f_18S_classified.tsv, results/step33f_composition.tsv.
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
DB=/srv/projects/db/SILVA/SSU/138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta
U=$W/results/step33f_18S_silva.tsv
if [ -e $U ]; then echo "REFUSING: $U exists"; exit 1; fi
hostname; date
$PY $W/scripts/step33f_collect.py || { echo "exit 1 (collect)"; exit 1; }
module load vsearch/2.21.1
module load workspace/scratch
vsearch --version 2>&1 | head -1
cp $DB $SCRATCH/silva.fa || { echo "ABORT: staging failed"; exit 1; }
vsearch --usearch_global $W/results/step33f_18S_ge500.fa --db $SCRATCH/silva.fa --id 0.70 --strand both \
  --maxaccepts 20 --maxrejects 100 --top_hits_only --notrunclabels --output_no_hits \
  --userout $SCRATCH/u.tsv --userfields query+target+id+alnlen+ql+qcov --threads $SLURM_CPUS_PER_TASK
rc=$?
echo "vsearch exit $rc"
if [ $rc -ne 0 ]; then exit 1; fi
cp $SCRATCH/u.tsv $U
$PY $W/scripts/step33f_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
