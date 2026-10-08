#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 8G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step14c_n88.log
# Step 14c n88: family sizes of the supported unknown reps at 88 metagenomes.
# No reclustering (n88 rule): uses results/unk_clusters_50_n88.tsv, the 50%
# identity, c 0.8 clustering of step14c restricted by step25b to members that
# are supported in the 88-metagenome catalog. Same member thresholds as step14c.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
F=$WPS/results/unk_clusters_50_n88.tsv
NU=$(awk -F'\t' 'NR>1 && $2=="1" && $6=="U"' $WPS/results/rep_classes_n88.tsv | wc -l)
NL=$(wc -l < $F)
echo "=== member lines: $NL, supported U reps n88: $NU ==="
[ "$NL" = "$NU" ] || { echo "GUARD FAIL: member lines differ from supported U count"; exit 1; }
awk -F'\t' -v N=$NU '{n[$1]++} END{
  for(r in n){t++; s=n[r];
    if(s>=3)a++; if(s>=25)b++; if(s>=50)c++; if(s>=100)d++;
    if(s>=3)pa+=s; if(s>=25)pb+=s; if(s>=50)pc+=s; if(s>=100)pd+=s}
  printf "families total\t%d\n", t;
  printf ">=3 members\t%d\tproteins %d (%.2f%%)\n", a, pa, 100*pa/N;
  printf ">=25 members\t%d\tproteins %d (%.2f%%)\n", b, pb, 100*pb/N;
  printf ">=50 members\t%d\tproteins %d (%.2f%%)\n", c, pc, 100*pc/N;
  printf ">=100 members\t%d\tproteins %d (%.2f%%)\n", d, pd, 100*pd/N
}' $F | tee $WPS/results/unk_family_sizes_n88.tsv
echo "DONE step14c_n88"
