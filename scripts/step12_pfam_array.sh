#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 24 --mem 32gb --time=3-00:00:00
#SBATCH --job-name=pfam
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/pfam.%A_%a.log
# Pfam domain assignment for the 95% representatives. Needed independently of
# eggNOG because the known-with-domain against known-without-domain split
# requires Pfam specifically. The database is staged to node-local storage and
# each task loops over every eighth chunk, matching the eggNOG job layout, since
# concurrent access to the shared copy stalled the eggNOG annotation phase.
module load hmmer/3.4
module load db-pfam
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
LOCAL=/scratch/lshad003_pfam_db
mkdir -p $LOCAL
for F in $PFAM_DB/Pfam-A.hmm*; do
  B=$(basename $F)
  [ -s "$LOCAL/$B" ] || cp "$F" "$LOCAL/$B"
done
echo "staged on $(hostname)"; ls -lh $LOCAL

IN=$WPS/catalog/db/LsPS_AA_95_rep__split
OUT=$WPS/results/pfam
mkdir -p $OUT
CPU=${SLURM_CPUS_ON_NODE:-24}
NTASKS=8
for (( N=SLURM_ARRAY_TASK_ID; N<=619; N+=NTASKS )); do
  F=$IN/LsPS.$N
  [ -s "$F" ] || continue
  [ -s "$OUT/LsPS.$N.domtblout" ] && { echo "skip $N"; continue; }
  echo "[$(date +%H:%M)] chunk $N"
  hmmscan --cut_ga --cpu $CPU \
    --domtblout $OUT/LsPS.$N.domtblout \
    --tblout $OUT/LsPS.$N.tblout \
    $LOCAL/Pfam-A.hmm $F > /dev/null
  echo "  hits: $(grep -vc '^#' $OUT/LsPS.$N.tblout 2>/dev/null || echo 0)"
done
echo "task $SLURM_ARRAY_TASK_ID done $(date)"
