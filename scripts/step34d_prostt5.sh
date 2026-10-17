#!/bin/bash
#SBATCH -J step34d_prostt5
#SBATCH -p stajichlab
#SBATCH --nodelist=r12
#SBATCH -c 32
#SBATCH --mem 64G
#SBATCH -t 2-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34d_prostt5.%j.log
# Step 34d: ProstT5 (CPU) + Foldseek structure screen of the 4,509 strict_final recurrent EU families, with 4,509
# length-matched strict_final K reps and 4,509 reversed decoys. Same binary, ProstT5 weights, databases and settings as
# step28h: createdb --prostt5-model, search vs PDB and CATH50 -e 1e-1, convertalis with the step28h fields.
# Inputs from scripts/step34d_inputs.py; summary scripts/step34d_summary.py. Output dir results/step34d_prostt5.
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
M=/bigdata/stajichlab/lshad003/software/foldseek_gpu/weights/prostt5
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
R=$W/results
O=$R/step34d_prostt5
if [ -e $O ]; then echo "REFUSING: $O exists"; exit 1; fi
hostname; date; $F version
$PY $W/scripts/step34d_inputs.py || { echo "exit 1 (inputs)"; exit 1; }
mkdir -p $O/tmp
SETS="$R/step34d_queries.faa:queries
$R/step34d_known_lenmatched.faa:known
$R/step34d_decoys_rev.faa:decoys"
for s in $SETS; do
  IFS=: read faa q <<< "$s"
  echo "== $q ($(grep -c '^>' $faa) sequences)"; date
  $F createdb $faa $O/$q --prostt5-model $M --threads $SLURM_CPUS_PER_TASK
  for t in PDB CATH50; do
    $F search $O/$q $DB/$t $O/${q}_res_$t $O/tmp -e 1e-1 --threads $SLURM_CPUS_PER_TASK
    $F convertalis $O/$q $DB/$t $O/${q}_res_$t $O/${q}_vs_$t.tsv \
      --format-output "query,target,evalue,bits,fident,alnlen,qlen,tlen,qcov,tcov"
    if [ ! -e $O/${q}_vs_$t.tsv ]; then echo "ABORT: no output for $q vs $t"; exit 1; fi
    echo "$q $t: hit rows $(wc -l < $O/${q}_vs_$t.tsv); queries with a hit $(cut -f1 $O/${q}_vs_$t.tsv | sort -u | wc -l)"
  done
done
date
$PY $W/scripts/step34d_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
