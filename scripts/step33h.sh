#!/bin/bash
#SBATCH -J step33h
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33h.%j.log
# Step 33h: scripts/step33h_classify.py (projection of EU reps on fly contigs onto dm6, GFF classes, repeat flags),
# DIAMOND 2.1.12 blastp --sensitive of results/step33h_queries.faa against the dm6 RefSeq proteome
# (catalog/fly/GCF_000001215.4_protein.faa.gz, 30,802 proteins) with --evalue 1e-3 --max-target-seqs 1, then
# scripts/step33h_summary.py. Outputs: results/step33h_rep_class.tsv, step33h_queries.faa, step33h_diamond_dm6.tsv,
# step33h_summary.tsv. Exits non-zero on any failure.
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
D=$W/results/step33h_diamond_dm6.tsv
if [ -e $D ]; then echo "REFUSING: $D exists"; exit 1; fi
hostname; date
$PY $W/scripts/step33h_classify.py || { echo "exit 1 (classify)"; exit 1; }
module load diamond/2.1.12
module load workspace/scratch
diamond version
diamond makedb --in $W/catalog/fly/GCF_000001215.4_protein.faa.gz -d $SCRATCH/dm6prot --threads $SLURM_CPUS_PER_TASK || { echo "exit 1 (makedb)"; exit 1; }
diamond blastp -q $W/results/step33h_queries.faa -d $SCRATCH/dm6prot -o $SCRATCH/dm.tsv --sensitive --evalue 1e-3 \
  --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --threads $SLURM_CPUS_PER_TASK || { echo "exit 1 (blastp)"; exit 1; }
cp $SCRATCH/dm.tsv $D
$PY $W/scripts/step33h_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
