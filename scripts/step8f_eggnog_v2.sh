#!/bin/bash
#SBATCH -p epyc -N 1 --cpus-per-task=24 --mem=120G --time=3-00:00:00
#SBATCH --job-name=eggnog2
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/eggnog2.%A_%a.log
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
SRC=/srv/projects/db/eggNOG/LATEST
LOCAL=/scratch/lshad003_eggnog_db
EMAPPER=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin/emapper.py
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin:$PATH

mkdir -p $LOCAL
for F in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
  [ -s "$LOCAL/$F" ] || cp "$SRC/$F" "$LOCAL/$F"
done
export EGGNOG_DATA_DIR=$LOCAL
echo "host: $(hostname)  db: $EGGNOG_DATA_DIR"
$EMAPPER --version

OUT=$WPS/results/eggnog2
mkdir -p $OUT $OUT/tmp_$SLURM_ARRAY_TASK_ID
NTASKS=8
for (( N=SLURM_ARRAY_TASK_ID; N<=619; N+=NTASKS )); do
  IN=$WPS/catalog/db/LsPS_AA_95_rep__split/LsPS.$N
  [ -s "$IN" ] || continue
  [ -s "$OUT/LsPS.$N.emapper.annotations" ] && \
    { [ "$(grep -vc '^#' $OUT/LsPS.$N.emapper.annotations)" -gt 0 ] && { echo "skip $N"; continue; }; }
  echo "[$(date +%H:%M)] chunk $N"
  $EMAPPER -i "$IN" -o LsPS.$N --output_dir "$OUT" \
    -m diamond --cpu 24 --itype proteins \
    --evalue 0.00001 --sensmode very-sensitive \
    --temp_dir $OUT/tmp_$SLURM_ARRAY_TASK_ID --override
  echo "  rows: $(grep -vc '^#' $OUT/LsPS.$N.emapper.annotations 2>/dev/null || echo FAIL)"
done
rm -rf $OUT/tmp_$SLURM_ARRAY_TASK_ID
echo "task $SLURM_ARRAY_TASK_ID done $(date)"
