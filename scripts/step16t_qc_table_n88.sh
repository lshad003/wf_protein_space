#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 4G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16t_n88.log
# Step 16t n88: mapping QC at 88 metagenomes. No remapping: the 89 table
# results/mapping_qc_89.tsv without UHM586.41010, written to results/mapping_qc_88.tsv.
R=/bigdata/stajichlab/lshad003/wf_protein_space/results
awk -F'\t' '$1!="UHM586.41010"' $R/mapping_qc_89.tsv > $R/mapping_qc_88.tsv
echo "qc rows: $(($(wc -l < $R/mapping_qc_88.tsv) - 1)), expect 88"
echo "=== mapped fraction by cohort ==="
awk -F'\t' 'NR>1{n[$2]++; s[$2]+=$5; if(mn[$2]==""||$5<mn[$2])mn[$2]=$5; if($5>mx[$2])mx[$2]=$5}
  END{for(c in n) printf "%s\tn=%d\tmean=%.4f\tmin=%s\tmax=%s\n",c,n[c],s[c]/n[c],mn[c],mx[c]}' \
  $R/mapping_qc_88.tsv | sort
echo "DONE step16t_n88"
