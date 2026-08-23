#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 24 --mem 48gb --time=3-00:00:00
#SBATCH --job-name=pfamS
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/pfamS.%A_%a.log
# Pfam domain assignment for the 95% representatives, using hmmsearch rather
# than hmmscan. hmmscan measured 78 minutes per 10,000-protein chunk, which
# exceeds the wall limit across 77 chunks per task; hmmsearch loads the profile
# database once and streams sequences past it. The two are equivalent searches
# but the table columns swap: in hmmsearch output the protein is the target
# (column 1) and the Pfam family is the query (column 3).
module load hmmer/3.4
module load db-pfam
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
LOCAL=/scratch/lshad003_pfam_db
mkdir -p $LOCAL
for F in $PFAM_DB/Pfam-A.hmm*; do
  B=$(basename $F); [ -s "$LOCAL/$B" ] || cp "$F" "$LOCAL/$B"
done
echo "staged on $(hostname) $(date)"

IN=$WPS/catalog/db/LsPS_AA_95_rep__split
OUT=$WPS/results/pfam_hs
mkdir -p $OUT
CPU=${SLURM_CPUS_ON_NODE:-24}
NTASKS=8
for (( N=SLURM_ARRAY_TASK_ID; N<=619; N+=NTASKS )); do
  F=$IN/LsPS.$N
  [ -s "$F" ] || continue
  [ -s "$OUT/LsPS.$N.done" ] && { echo "skip $N"; continue; }
  echo "[$(date +%H:%M)] chunk $N"
  hmmsearch --cut_ga --cpu $CPU \
    --domtblout $OUT/LsPS.$N.domtblout \
    --tblout $OUT/LsPS.$N.tblout \
    $LOCAL/Pfam-A.hmm $F > /dev/null
  touch $OUT/LsPS.$N.done
  echo "  proteins with a domain: $(grep -v '^#' $OUT/LsPS.$N.tblout | awk '{print $1}' | sort -u | wc -l)"
done
echo "task $SLURM_ARRAY_TASK_ID done $(date)"
