#!/bin/bash
#SBATCH -p stajichlab -c 16 --mem 32G --time 8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step14b.log
# Step 14b: AntiFam screen of the 586,215 supported unknown reps.
# AntiFam flags spurious ORFs, pseudogenes, false translations.
# Same control the Pavlopoulos 2023 paper used (they found 43 sequences).
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
Q=$WPS/results/supported_U_reps.faa
AF=$WPS/db/antifam

N=$(grep -c '>' $Q)
[ "$N" = "586215" ] || { echo "GUARD FAIL: query $N need 586215"; exit 1; }

mkdir -p $AF
if [ ! -s $AF/AntiFam.hmm ]; then
  cd $AF
  curl -L -o Antifam.tar.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz || exit 1
  tar xzf Antifam.tar.gz || exit 1
  ls -l
fi
[ -s $AF/AntiFam.hmm ] || { echo "AntiFam.hmm not found after download:"; ls $AF; exit 1; }

module load hmmer/3.4 2>/dev/null || module load hmmer/3.3.2
hmmpress -f $AF/AntiFam.hmm > /dev/null 2>&1
echo "hmmer: $(hmmsearch -h | head -2 | tail -1)"
echo "AntiFam models: $(grep -c '^NAME' $AF/AntiFam.hmm)"

hmmsearch --cut_ga --cpu 16 \
  --tblout $WPS/results/antifam_hits.tblout \
  -o /dev/null \
  $AF/AntiFam.hmm $Q

HITS=$(grep -v '^#' $WPS/results/antifam_hits.tblout | awk '{print $1}' | sort -u | wc -l)
echo "=== spurious sequences flagged: $HITS of 586215 ($(awk -v h=$HITS 'BEGIN{printf "%.4f", 100*h/586215}')%) ==="
grep -v '^#' $WPS/results/antifam_hits.tblout | awk '{print $3}' | sort | uniq -c | sort -rn | head -10
echo "DONE step14b"
