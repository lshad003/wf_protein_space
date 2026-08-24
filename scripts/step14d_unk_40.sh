#!/bin/bash
#SBATCH -p epyc -c 24 --mem 120G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step14d.log
# Step 14d: cluster unknown reps at 40% identity, sensitivity -s 7.5,
# to test whether the non-clustering at 50% is a threshold effect.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
D=$WPS/catalog/db/UNK
mkdir -p $WPS/tmp/unk40
module load mmseqs2/13-45111
mmseqs cluster $D/UNK_AA $D/UNK_40_cluster $WPS/tmp/unk40 \
  --min-seq-id 0.4 -c 0.8 --cov-mode 0 -s 7.5 --threads 24
mmseqs createtsv $D/UNK_AA $D/UNK_AA $D/UNK_40_cluster $WPS/results/unk_clusters_40.tsv --threads 24
awk -F'\t' '{n[$1]++} END{
  for(r in n){t++; s=n[r]; if(s>=3)a++; if(s>=25)b++; if(s>=100)d++}
  printf "40%% identity: families %d | >=3 %d | >=25 %d | >=100 %d\n", t,a,b,d
}' $WPS/results/unk_clusters_40.tsv
echo "DONE step14d"
