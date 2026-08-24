#!/bin/bash
# Step 21a: confirm gene IDs encode contig and gene number, and that
# neighbours of unknown genes exist in the catalog with annotations.
W=/bigdata/stajichlab/lshad003/wf_protein_space
export LC_ALL=C

echo "=== 1. gene ID format, a few unknown reps ==="
head -3 $W/results/supported_U_reps.faa | grep '>'
awk -F'\t' 'NR>1 && $2=="EU"{print $1; c++} c==3{exit}' $W/results/fourway_classes.tsv

echo ""
echo "=== 2. do neighbours appear in clusters_95.tsv (any member, not just reps)? ==="
G=$(awk -F'\t' 'NR>1 && $2=="EU"{print $1; exit}' $W/results/fourway_classes.tsv)
echo "gene: $G"
STEM=${G%%__*}; REST=${G#*__}; NUM=${REST##*_}; CONTIG=${REST%_*}
echo "stem=$STEM contig=$CONTIG num=$NUM"
for N in $((NUM-1)) $NUM $((NUM+1)); do
  echo -n "  ${STEM}__${CONTIG}_${N}: "
  grep -c -F "${STEM}__${CONTIG}_${N}" $W/results/clusters_95.tsv 2>/dev/null || echo 0
done

echo ""
echo "=== 3. unknown families with >=3 members, sizes ==="
awk -F'\t' '{n[$1]++} END{for(r in n) if(n[r]>=3) print n[r]}' $W/results/unk_clusters_50.tsv \
| sort -rn | head -5
echo "families >=3: $(awk -F'\t' '{n[$1]++} END{c=0; for(r in n) if(n[r]>=3)c++; print c}' $W/results/unk_clusters_50.tsv)"
