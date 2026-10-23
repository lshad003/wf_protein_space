#!/bin/bash
#SBATCH -J step29b_esmfold
#SBATCH -p epyc
#SBATCH -c 32
#SBATCH --mem 128G
#SBATCH -t 7-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step29b_esmfold.%j.log
# Step 29b: ESMFold (esm/2.0.0, esmfold_3B_v1) on results/esmfold_input.faa (218 HIT + 200 NOHIT + 200 NOHIT_LONG, step29a).
# One PDB per sequence in results/esmfold_pdb; pLDDT in the B-factor column. --chunk-size 64 so the
# 1,077 aa query fits in memory.
# CPU only on epyc: GPU attempts never started (gpu job 29931830 est. 2026-10-25; preempt_gpu A100 jobs
# 29931855, 29931930 held for higher-priority partitions), all cancelled.
# Torch threads limited to the 32 allocated cores. Refuses to run if the output directory exists.
W=/bigdata/stajichlab/lshad003/wf_protein_space
OUT=$W/results/esmfold_pdb
if [ -e $OUT ]; then echo "REFUSING: $OUT exists"; exit 1; fi
export TORCH_HOME=/srv/projects/db/esm
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
hostname; date
/opt/linux/rocky/8.x/x86_64/pkgs/esm/2.0.0/env/bin/esm-fold -i $W/results/esmfold_input.faa -o $OUT --chunk-size 64 --cpu-only
rc=$?
echo "exit $rc"; date; ls $OUT | grep -c "\.pdb$"
# 2026-10-23: exit with esm-fold's status so afterok dependencies see failures (CLAUDE.md rule)
exit $rc
