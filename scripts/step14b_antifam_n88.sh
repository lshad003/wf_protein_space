#!/bin/bash
#SBATCH -p stajichlab -c 16 --mem 32G --time 8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step14b_n88.log
# Step 14b n88: AntiFam screen of the supported unknown reps of the
# 88-metagenome catalog. Query is regenerated as results/supported_U_reps_n88.faa:
# the records of supported_U_reps.faa whose IDs are U and supported in
# results/rep_classes_n88.tsv. Same AntiFam copy (db/antifam) and --cut_ga.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
Q=$WPS/results/supported_U_reps_n88.faa
AF=$WPS/db/antifam
[ -s $AF/AntiFam.hmm ] || { echo "GUARD FAIL: $AF/AntiFam.hmm missing"; exit 1; }

awk -F'\t' 'NR>1 && $2=="1" && $6=="U"{print $1}' $WPS/results/rep_classes_n88.tsv > $WPS/results/supported_U_ids_n88.txt
NU=$(wc -l < $WPS/results/supported_U_ids_n88.txt)
awk 'NR==FNR{k[$1]=1; next} /^>/{id=substr($1,2); p=(id in k)} p' \
  $WPS/results/supported_U_ids_n88.txt $WPS/results/supported_U_reps.faa > $Q
N=$(grep -c '>' $Q)
echo "supported U ids n88: $NU  records written: $N"
[ "$N" = "$NU" ] && [ "$N" = "585347" ] || { echo "GUARD FAIL: query $N, ids $NU, need 585347"; exit 1; }

module load hmmer/3.4 2>/dev/null || module load hmmer/3.3.2
echo "hmmer: $(hmmsearch -h | head -2 | tail -1)"
echo "AntiFam models: $(grep -c '^NAME' $AF/AntiFam.hmm)  version: $(cat $AF/version)"

hmmsearch --cut_ga --cpu 16 \
  --tblout $WPS/results/antifam_hits_n88.tblout \
  -o /dev/null \
  $AF/AntiFam.hmm $Q

HITS=$(grep -v '^#' $WPS/results/antifam_hits_n88.tblout | awk '{print $1}' | sort -u | wc -l)
echo "=== spurious sequences flagged: $HITS of $N ($(awk -v h=$HITS -v n=$N 'BEGIN{printf "%.4f", 100*h/n}')%) ==="
grep -v '^#' $WPS/results/antifam_hits_n88.tblout | awk '{print $3}' | sort | uniq -c | sort -rn | head -10
OLDSUB=$(grep -v '^#' $WPS/results/antifam_hits.tblout | awk '{print $1}' | sort -u | grep -c -x -F -f $WPS/results/supported_U_ids_n88.txt)
echo "89-run hits whose rep is in the n88 query: $OLDSUB (expect equal to $HITS)"
echo "DONE step14b_n88"
