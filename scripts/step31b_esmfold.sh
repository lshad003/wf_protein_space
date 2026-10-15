#!/bin/bash
#SBATCH -J step31b_esmfold
#SBATCH -p epyc
#SBATCH -c 32
#SBATCH --mem 128G
#SBATCH -t 2-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step31b_esmfold.%j.log
# Step 31b: ESMFold (esm/2.0.0, esmfold_3B_v1) of the 73 defence-context single-variant reps not folded in Step 29
# (results/def84_fold.faa, step31a). Same call as step29b: --chunk-size 64, --cpu-only, TORCH_HOME, threads limited
# to the allocated cores. One PDB per sequence in results/esmfold_def84_pdb; pLDDT in the B-factor column.
# Refuses to run if the output directory exists.
W=/bigdata/stajichlab/lshad003/wf_protein_space
OUT=$W/results/esmfold_def84_pdb
if [ -e $OUT ]; then echo "REFUSING: $OUT exists"; exit 1; fi
export TORCH_HOME=/srv/projects/db/esm
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
hostname; date
/opt/linux/rocky/8.x/x86_64/pkgs/esm/2.0.0/env/bin/esm-fold -i $W/results/def84_fold.faa -o $OUT --chunk-size 64 --cpu-only
echo "exit $?"; date; ls $OUT | grep -c "\.pdb$"
