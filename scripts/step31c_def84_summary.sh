#!/bin/bash
#SBATCH -J step31c_def84
#SBATCH -p epyc
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 4:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step31c_def84_summary.%j.log
# Step 31c: local Foldseek of the 84 defence-context models and the per-model summary.
# Queries: results/esmfold_def84_pdb (step31b, 73) + results/def84_step29_pdb (Step 29 copies, 11) collected in
# results/foldseek_def84/; aborts unless 84.
# Binary and databases as step29e: Foldseek commit 463739e0014a1549a527de589102cde98f802f37 on CPU; PDB and CATH50
# from /bigdata/stajichlab/lshad003/software/foldseek_dbs. easy-search --exhaustive-search 1 -e 10, step29e fields
# plus theader. Summary in scripts/step31c_def84_summary.py.
# Outputs (refuses to overwrite): results/foldseek_def84/, results/step31_foldseek_{PDB,CATH50}.m8,
# results/step31_def84_summary.tsv
W=/bigdata/stajichlab/lshad003/wf_protein_space
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
Q=$W/results/foldseek_def84
for p in $Q $W/results/step31_def84_summary.tsv $W/results/step31_foldseek_PDB.m8 $W/results/step31_foldseek_CATH50.m8; do
  if [ -e $p ]; then echo "REFUSING: $p exists"; exit 1; fi
done
hostname; date; $F version
mkdir -p $Q
cp $W/results/esmfold_def84_pdb/*.pdb $W/results/def84_step29_pdb/*.pdb $Q/
n=$(ls $Q | grep -c '\.pdb$')
echo "queries: $n"
if [ "$n" -ne 84 ]; then echo "ABORT: expected 84 query PDBs, got $n"; exit 1; fi
module load workspace/scratch
T=$SCRATCH/step31c
mkdir -p $T
for db in PDB CATH50; do
  echo "== $db"; date
  $F easy-search $Q $DB/$db $T/$db.m8 $T/tmp_$db --exhaustive-search 1 -e 10 --threads $SLURM_CPUS_PER_TASK \
    --format-output query,target,theader,fident,alnlen,qlen,tlen,evalue,bits,prob,alntmscore,qtmscore,ttmscore,lddt,qcov,tcov
  echo "exit $?"
  if [ ! -e $T/$db.m8 ]; then echo "ABORT: no output for $db"; rm -rf $T; exit 1; fi
  cp $T/$db.m8 $W/results/step31_foldseek_$db.m8
  echo "$db: rows $(wc -l < $W/results/step31_foldseek_$db.m8), queries with a hit $(cut -f1 $W/results/step31_foldseek_$db.m8 | sort -u | wc -l)"
done
rm -rf $T
$PY $W/scripts/step31c_def84_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
