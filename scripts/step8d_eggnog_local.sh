#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 64gb --time=3-00:00:00
#SBATCH --job-name=egg_local
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/egglocal.%A_%a.log
module load eggnog-mapper/2.1.9
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
SRC=/bigdata/operations/pkgadmin/srv/projects/db/eggNOG/5.0.2
INDIR=$WPS/catalog/db/LsPS_AA_95_rep__split
OUTDIR=$WPS/results/eggnog
mkdir -p $OUTDIR

LOCAL=/scratch/lshad003_eggnog_db
mkdir -p $LOCAL
echo "staging eggNOG db to $LOCAL on $(hostname) at $(date)"
for F in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
  if [ ! -s "$LOCAL/$F" ] || [ "$(stat -c %s $SRC/$F)" != "$(stat -c %s $LOCAL/$F 2>/dev/null || echo 0)" ]; then
    echo "  copying $F"; cp "$SRC/$F" "$LOCAL/$F"
  else
    echo "  $F already staged"
  fi
done
echo "staging done at $(date)"
df -h /scratch | tail -1

# each array task takes every Nth chunk
NTASKS=8
T=$SLURM_ARRAY_TASK_ID
for (( N=T; N<=619; N+=NTASKS )); do
  FASTA=$INDIR/LsPS.$N
  [ -s "$FASTA" ] || continue
  [ -s "$OUTDIR/LsPS.$N.emapper.annotations" ] && { echo "skip $N (done)"; continue; }
  HITS=$OUTDIR/LsPS.$N.emapper.hits
  if [ -s "$HITS" ]; then
    echo "[$(date +%H:%M)] chunk $N: annotating from existing hits"
    emapper.py -m no_search --annotate_hits_table "$HITS" \
      -o LsPS.$N --output_dir "$OUTDIR" --cpu 32 --data_dir "$LOCAL" --override
  else
    echo "[$(date +%H:%M)] chunk $N: full run"
    emapper.py -i "$FASTA" -o LsPS.$N --output_dir "$OUTDIR" \
      --cpu 32 -m diamond --itype proteins \
      --evalue 0.00001 --sensmode very-sensitive --data_dir "$LOCAL"
  fi
done
echo "task $T finished at $(date)"
