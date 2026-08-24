#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 4G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16t.log
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
R=$WPS/results
M=$WPS/metadata

awk '{print $1"\tWF22"}' $M/wf22_stems_44.txt  > $R/cohort_map.tsv
awk '{print $1"\tWF23"}' $M/wf23_stems.txt    >> $R/cohort_map.tsv
awk '{print $1"\tWF24"}' $M/wf24_stems.txt    >> $R/cohort_map.tsv
echo "cohort map rows: $(wc -l < $R/cohort_map.tsv), expect 89"

BAD=0
printf "stem\tcohort\trecords\tmapped_primary\tfraction\n" > $R/mapping_qc_89.tsv
while IFS=$'\t' read -r S C; do
  F=$R/counts/$S.counts.tsv.gz
  U=$R/mapping_summary/$S.summary.tsv
  if [ ! -s "$F" ] || ! gzip -t "$F" 2>/dev/null; then echo "BAD counts: $S"; BAD=$((BAD+1)); continue; fi
  if [ "$(awk 'NF==4' $U 2>/dev/null | wc -l)" != "1" ]; then echo "BAD summary: $S"; BAD=$((BAD+1)); continue; fi
  awk -F'\t' -v c="$C" '{printf "%s\t%s\t%s\t%s\t%s\n",$1,c,$2,$3,$4}' $U >> $R/mapping_qc_89.tsv
done < $R/cohort_map.tsv
echo "bad files: $BAD, expect 0"
echo "qc rows: $(($(wc -l < $R/mapping_qc_89.tsv) - 1)), expect 89"

echo "=== mapped fraction by cohort ==="
awk -F'\t' 'NR>1{n[$2]++; s[$2]+=$5; if(mn[$2]==""||$5<mn[$2])mn[$2]=$5; if($5>mx[$2])mx[$2]=$5}
  END{for(c in n) printf "%s\tn=%d\tmean=%.4f\tmin=%s\tmax=%s\n",c,n[c],s[c]/n[c],mn[c],mx[c]}' \
  $R/mapping_qc_89.tsv | sort
echo "DONE step16t"
