#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 4 --mem 8gb --time=2-00:00:00
#SBATCH --job-name=pfam
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/pfam.%A_%a.log
module load hmmer/3.4
module load db-pfam
module load workspace/scratch
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
CPU=${SLURM_CPUS_ON_NODE:-4}
N=$SLURM_ARRAY_TASK_ID
IN=$WPS/catalog/db/LsPS_AA_95_rep__split/LsPS.$N
OUT=$WPS/results/pfam
mkdir -p $OUT
[ -s "$IN" ] || { echo "no chunk $N"; exit 0; }
[ -s "$OUT/LsPS.$N.domtblout" ] && { echo "done already"; exit 0; }
rsync -a $PFAM_DB/Pfam-A.hmm* $SCRATCH/
hmmscan --cut_ga --cpu $CPU \
  --domtblout $OUT/LsPS.$N.domtblout \
  --tblout $OUT/LsPS.$N.tblout \
  $SCRATCH/Pfam-A.hmm $IN > /dev/null
echo "finished chunk $N"
