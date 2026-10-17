#!/bin/bash
#SBATCH -J step34d2_kfix
#SBATCH -p stajichlab
#SBATCH --nodelist=r12
#SBATCH -c 32
#SBATCH --mem 64G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34d2_kfix.%j.log
# Step 34d2: ProstT5/Foldseek on the 225 replacement K reps only, same settings as step34d_prostt5.sh, then merge.
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
M=/bigdata/stajichlab/lshad003/software/foldseek_gpu/weights/prostt5
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
O=$W/results/step34d2_prostt5
if [ -e $O ]; then echo "REFUSING: $O exists"; exit 1; fi
hostname; date; $F version
$PY $W/scripts/step34d2_kfix_inputs.py || { echo "exit 1 (inputs)"; exit 1; }
mkdir -p $O/tmp
$F createdb $W/results/step34d2_k_replacements.faa $O/krepl --prostt5-model $M --threads $SLURM_CPUS_PER_TASK
for t in PDB CATH50; do
  $F search $O/krepl $DB/$t $O/krepl_res_$t $O/tmp -e 1e-1 --threads $SLURM_CPUS_PER_TASK
  $F convertalis $O/krepl $DB/$t $O/krepl_res_$t $O/krepl_vs_$t.tsv \
    --format-output "query,target,evalue,bits,fident,alnlen,qlen,tlen,qcov,tcov"
  if [ ! -e $O/krepl_vs_$t.tsv ]; then echo "ABORT: no output for $t"; exit 1; fi
  echo "krepl $t: hit rows $(wc -l < $O/krepl_vs_$t.tsv); queries with a hit $(cut -f1 $O/krepl_vs_$t.tsv | sort -u | wc -l)"
done
echo "createdb entries: $(wc -l < $O/krepl.index)"
$PY $W/scripts/step34d2_kfix_merge.py
rc=$?
echo "exit $rc"; date
exit $rc
