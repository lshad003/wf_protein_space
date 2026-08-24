#!/bin/bash
# Step 20b: check the LCA file format before writing the join.
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
S=UHM102.10840
echo "=== header and first lines ==="
head -3 $FEC/results_scaffold_classify_mmseqs/$S/${S}_uniref50_lca.tsv
echo ""
echo "=== columns ==="
head -1 $FEC/results_scaffold_classify_mmseqs/$S/${S}_uniref50_lca.tsv | awk -F'\t' '{print NF" columns"}'
echo ""
echo "=== how many stems have this file ==="
W=/bigdata/stajichlab/lshad003/wf_protein_space
N=0
for X in $(cat $W/metadata/wf22_stems_44.txt $W/metadata/wf23_stems.txt $W/metadata/wf24_stems.txt); do
  [ -s $FEC/results_scaffold_classify_mmseqs/$X/${X}_uniref50_lca.tsv ] && N=$((N+1))
done
echo "$N of 89"
