#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 8G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step13h_n88.log
# Step 13h n88: GU/EU split under the step13h criteria, with hits restricted to
# the supported U reps of the 88-metagenome catalog (results/supported_U_ids_n88.txt,
# written by step14b_n88) and the denominator set to that count. Same search output.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
H=$WPS/results/step13_nrclust20260128_hits.tsv
IDS=$WPS/results/supported_U_ids_n88.txt
export LC_ALL=C
T=$(wc -l < $IDS)
[ "$T" = "585347" ] || { echo "GUARD FAIL: $IDS has $T, need 585347"; exit 1; }
awk -F'\t' 'NR==FNR{k[$1]=1; next} ($1 in k)' $IDS $H > $WPS/tmp/step13h_n88_hits.tsv
echo "hit lines: 89 $(wc -l < $H), n88 $(wc -l < $WPS/tmp/step13h_n88_hits.tsv); supported U n88: $T"
echo ""
echo "criterion                          GU(hit)   EU(no hit)   EU%"
for C in "1e-5:0:0" "1e-5:50:0" "1e-10:50:50" "1e-5:70:70" "1e-20:50:50"; do
  E=$(echo $C | cut -d: -f1); QC=$(echo $C | cut -d: -f2); SC=$(echo $C | cut -d: -f3)
  N=$(awk -F'\t' -v e=$E -v q=$QC -v s=$SC '$6<=e && $8>=q && $9>=s {print $1}' $WPS/tmp/step13h_n88_hits.tsv | sort -u | wc -l)
  printf "E<=%-7s qcov>=%-3s scov>=%-3s   %-9d %-12d %.2f%%\n" $E $QC $SC $N $((T-N)) $(awk -v a=$((T-N)) -v t=$T 'BEGIN{printf "%.2f", 100*a/t}')
done
rm $WPS/tmp/step13h_n88_hits.tsv
echo "DONE step13h_n88"
