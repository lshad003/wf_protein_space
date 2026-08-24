#!/bin/bash
# Step 13h: GU/EU split under several stated criteria. Login-cheap, awk only.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
H=$WPS/results/step13_nrclust20260128_hits.tsv
export LC_ALL=C
T=586215
echo "hit lines: $(wc -l < $H)"
echo ""
echo "criterion                          GU(hit)   EU(no hit)   EU%"
for C in "1e-5:0:0" "1e-5:50:0" "1e-10:50:50" "1e-5:70:70" "1e-20:50:50"; do
  E=$(echo $C | cut -d: -f1); QC=$(echo $C | cut -d: -f2); SC=$(echo $C | cut -d: -f3)
  N=$(awk -F'\t' -v e=$E -v q=$QC -v s=$SC '$6<=e && $8>=q && $9>=s {print $1}' $H | sort -u | wc -l)
  printf "E<=%-7s qcov>=%-3s scov>=%-3s   %-9d %-12d %.2f%%\n" $E $QC $SC $N $((T-N)) $(awk -v a=$((T-N)) -v t=$T 'BEGIN{printf "%.2f", 100*a/t}')
done
echo ""
echo "=== top subject titles among hits (what the GU genes resemble) ==="
awk -F'\t' '$6<=1e-5 {print $10}' $H | sed 's/ \[.*//' | sort | uniq -c | sort -rn | head -12
