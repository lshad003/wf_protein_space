#!/bin/bash
#SBATCH -J step34j_esmfold
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem 96G
#SBATCH -t 3-00:00:00
#SBATCH --array=0-7
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_b1_esmfold.%A_%a.log
# Step 34j B1: ESMFold (esm/2.0.0 module, esmfold_3B_v1; Step 29b setup) on results/step34j/esmfold_chunks/C.<task>.faa
# (8,759 sequences: 4,250 member-validated EU family reps + 4,509 known controls; scripts/step34j_b1_esmfold_inputs.py).
# GPU instead of the Step 29b CPU run; --chunk-size 64 as Step 29b so the 2,839 aa query fits.
# One PDB per sequence in results/step34j/esmfold_pdb (pLDDT in the B-factor column); mean pLDDT and pTM are in this
# log ("Predicted structure for ... pLDDT x, pTM y"). Refuses to run if any PDB of this chunk already exists.
W=/bigdata/stajichlab/lshad003/wf_protein_space
I=$W/results/step34j/esmfold_chunks/C.$SLURM_ARRAY_TASK_ID.faa
OUT=$W/results/step34j/esmfold_pdb
if [ ! -s $I ]; then echo "ABORT: no input $I"; exit 1; fi
mkdir -p $OUT
for id in $(grep '^>' $I | cut -c2-); do
  if [ -e $OUT/$id.pdb ]; then echo "REFUSING: $OUT/$id.pdb exists"; exit 1; fi
done
export TORCH_HOME=/srv/projects/db/esm
hostname; date; nvidia-smi --query-gpu=name,memory.total --format=csv
echo "sequences $(grep -c '^>' $I)"
/opt/linux/rocky/8.x/x86_64/pkgs/esm/2.0.0/env/bin/esm-fold -i $I -o $OUT --chunk-size 64
rc=$?
echo "exit $rc"; date
n=0; for id in $(grep '^>' $I | cut -c2-); do [ -s $OUT/$id.pdb ] && n=$((n+1)); done
echo "models written for this chunk: $n of $(grep -c '^>' $I)"
# 2026-10-23: exit with esm-fold's status so afterok dependencies see failures (CLAUDE.md rule)
exit $rc
