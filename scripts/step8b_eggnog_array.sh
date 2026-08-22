#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -c 4 --mem 16gb --time=2-00:00:00
#SBATCH --job-name=eggnog
#SBATCH --output=/bigdata/stajichlab/lshad003/wf_protein_space/logs/eggnog.%A_%a.log
module load eggnog-mapper/2.1.9
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
INDIR=$WPS/catalog/db/LsPS_AA_95_rep__split
OUTDIR=$WPS/results/eggnog
mkdir -p $OUTDIR

N=$((SLURM_ARRAY_TASK_ID))
FASTA=$INDIR/LsPS.$N
[ -s "$FASTA" ] || { echo "no chunk $N"; exit 0; }
PREFIX=LsPS.$N

if [ -s "$OUTDIR/${PREFIX}.emapper.annotations" ]; then
  echo "already done: $PREFIX"; exit 0
fi

emapper.py -i "$FASTA" -o "$PREFIX" --output_dir "$OUTDIR" \
  --cpu 4 -m diamond --itype proteins \
  --evalue 0.00001 --sensmode very-sensitive
echo "finished $FASTA"
