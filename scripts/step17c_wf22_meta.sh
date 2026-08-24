#!/bin/bash
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
OUT=$WPS/metadata/wf22_design.tsv
export LC_ALL=C

awk -F'\t' 'NR>1 && $3!=""{split($1,a,"."); print $1"\t"a[1]"\t"$3"\t"$5"\t"$2}' \
  $V1/last_wood_frog_2022_sample_metadata.tsv | sort > /tmp/m.tsv
sort $WPS/metadata/wf22_stems_44.txt > /tmp/ids.txt
printf "stem\tanimal\ttreatment\tegg_mass\tmonth\n" > $OUT
join -t$'\t' /tmp/ids.txt /tmp/m.tsv >> $OUT

echo "rows: $(($(wc -l < $OUT) - 1)), expect 44"
echo ""
echo "=== samples per animal ==="
awk -F'\t' 'NR>1{print $2}' $OUT | sort | uniq -c | awk '{n[$1]++} END{for(k in n) print k" samples: "n[k]" animals"}' | sort
echo ""
echo "=== treatment x egg_mass ==="
awk -F'\t' 'NR>1{print $3"\t"$4}' $OUT | sort | uniq -c
echo ""
echo "=== month ==="
awk -F'\t' 'NR>1{print $5}' $OUT | sort | uniq -c
