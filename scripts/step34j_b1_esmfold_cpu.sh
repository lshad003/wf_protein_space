#!/bin/bash
#SBATCH -J step34j_esmfold_cpu
#SBATCH -p stajichlab
#SBATCH -c 28
#SBATCH --mem 110G
#SBATCH -t 10-00:00:00
#SBATCH --array=0-23
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_b1_esmfold_cpu.%A_%a.log
# Step 34j B1 (CPU): ESMFold (esm/2.0.0 module, esmfold_3B_v1; exactly the Step 29b setup, --cpu-only, --chunk-size 64)
# on results/step34j/esmfold_chunks_cpu/C.<task>.faa (8,759 sequences: 4,250 member-validated EU family reps + 4,509
# known controls). Replaces the GPU array (scripts/step34j_b1_esmfold.sh, job 30049857, cancelled unstarted).
# One PDB per sequence in results/step34j/esmfold_pdb; mean pLDDT and pTM in this log. Refuses if any PDB of this
# chunk already exists.
W=/bigdata/stajichlab/lshad003/wf_protein_space
I=$W/results/step34j/esmfold_chunks_cpu/C.$SLURM_ARRAY_TASK_ID.faa
OUT=$W/results/step34j/esmfold_pdb
if [ ! -s $I ]; then echo "ABORT: no input $I"; exit 1; fi
mkdir -p $OUT
for id in $(grep '^>' $I | cut -c2-); do
  if [ -e $OUT/$id.pdb ]; then echo "REFUSING: $OUT/$id.pdb exists"; exit 1; fi
done
export TORCH_HOME=/srv/projects/db/esm
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
hostname; date; echo "sequences $(grep -c '^>' $I)"
/opt/linux/rocky/8.x/x86_64/pkgs/esm/2.0.0/env/bin/esm-fold -i $I -o $OUT --chunk-size 64 --cpu-only
echo "exit $?"; date
n=0; for id in $(grep '^>' $I | cut -c2-); do [ -s $OUT/$id.pdb ] && n=$((n+1)); done
echo "models written for this chunk: $n of $(grep -c '^>' $I)"
