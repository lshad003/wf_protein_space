#!/bin/bash
# Step 16c: locate nucleotide CDS and clean reads by direct listing.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal

echo "=== 1. step1c variables, where the 14 CDS files went ==="
grep -nE '^OUT=|^IN=|^DIR=|BASE=' $BASE/scripts/step1c_prodigal_array.sh

echo ""
echo "=== 2. step1d rename transform, first 30 lines ==="
head -30 $BASE/scripts/step1d_rename_45.sh

echo ""
echo "=== 3. one per-stem dir per cohort under Fecal/results ==="
S22=$(head -1 $BASE/metadata/wf22_stems_44.txt 2>/dev/null || head -1 $BASE/metadata/wf22_stems.txt)
S23=$(head -1 $BASE/metadata/wf23_stems.txt)
S24=$(head -1 $BASE/metadata/wf24_stems.txt)
for S in $S22 $S23 $S24; do
  echo "--- $FEC/results/$S ---"
  ls $FEC/results/$S 2>/dev/null | head -25
done

echo ""
echo "=== 4. top-level dirs of the project, candidate CDS locations ==="
ls -d $BASE/*/ 2>/dev/null

echo ""
echo "=== 5. WF24 treatment codes, resolves 7x5=35 vs 36 ==="
awk -F'\t' '$3=="WF24"{print $5}' $BASE/metadata/wf23_wf24_sequenced_metadata.tsv | sort | uniq -c
