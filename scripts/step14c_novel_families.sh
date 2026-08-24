#!/bin/bash
#SBATCH -p epyc -c 24 --mem 120G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step14c.log
# Step 14c: cluster the 586,215 supported unknown reps into families at 50%
# identity, then count families at member thresholds comparable to
# Pavlopoulos 2023 (>=3, >=25, >=50, >=100 members).
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
Q=$WPS/results/supported_U_reps.faa
D=$WPS/catalog/db/UNK
mkdir -p $D $WPS/tmp/unk

N=$(grep -c '>' $Q)
[ "$N" = "586215" ] || { echo "GUARD FAIL: query $N need 586215"; exit 1; }

module load mmseqs2/13-45111
mmseqs version

mmseqs createdb $Q $D/UNK_AA
mmseqs cluster $D/UNK_AA $D/UNK_50_cluster $WPS/tmp/unk \
  --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads 24
mmseqs createtsv $D/UNK_AA $D/UNK_AA $D/UNK_50_cluster $WPS/results/unk_clusters_50.tsv --threads 24

echo "=== member lines: $(wc -l < $WPS/results/unk_clusters_50.tsv), expect 586215 ==="
awk -F'\t' '{n[$1]++} END{
  for(r in n){t++; s=n[r];
    if(s>=3)a++; if(s>=25)b++; if(s>=50)c++; if(s>=100)d++;
    if(s>=3)pa+=s; if(s>=25)pb+=s; if(s>=50)pc+=s; if(s>=100)pd+=s}
  printf "families total\t%d\n", t;
  printf ">=3 members\t%d\tproteins %d (%.2f%%)\n", a, pa, 100*pa/586215;
  printf ">=25 members\t%d\tproteins %d (%.2f%%)\n", b, pb, 100*pb/586215;
  printf ">=50 members\t%d\tproteins %d (%.2f%%)\n", c, pc, 100*pc/586215;
  printf ">=100 members\t%d\tproteins %d (%.2f%%)\n", d, pd, 100*pd/586215
}' $WPS/results/unk_clusters_50.tsv | tee $WPS/results/unk_family_sizes.tsv
echo "DONE step14c"
