#!/bin/bash
# Step 16i: fastq locations via manifests and v1 root scripts; step16h check.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1

echo "=== 1. step16h state and log tail ==="
sacct -j 27713792 -n -X --format=State,Elapsed,MaxRSS
tail -6 $BASE/logs/step16h.log 2>/dev/null

echo ""
echo "=== 2. manifests at Fecal top level ==="
for F in read_manifest.csv samples.csv names_w_path.txt missing_reads.txt; do
  echo "--- $F ($(wc -l < $FEC/$F 2>/dev/null) lines) ---"
  head -3 $FEC/$F 2>/dev/null
done

echo ""
echo "=== 3. read dirs, one stem ==="
echo "--- all_input_filtered ---"
ls $FEC/all_input_filtered 2>/dev/null | head -5
ls $FEC/all_input_filtered/UHM102.10840* 2>/dev/null
echo "--- input ---"
ls $FEC/input 2>/dev/null | head -5
ls $FEC/input/UHM102.10840* 2>/dev/null

echo ""
echo "=== 4. v1 root mapping scripts, read path lines ==="
for F in map_basidiobolus_minimap2_array.sh map_basidiobolus_diamond_array.sh; do
  echo "--- $F ---"
  grep -nE 'fq|fastq|R1|R2|input|filtered' $V1/$F 2>/dev/null | head -8
done

echo ""
echo "=== 5. presence across all 89 stems ==="
MISS=0
for S in $(cat $BASE/metadata/wf22_stems_44.txt $BASE/metadata/wf23_stems.txt $BASE/metadata/wf24_stems.txt); do
  ok=0
  ls $FEC/all_input_filtered/${S}* > /dev/null 2>&1 && ok=1
  [ $ok -eq 0 ] && ls $FEC/input/${S}* > /dev/null 2>&1 && ok=2
  [ $ok -eq 0 ] && { MISS=$((MISS+1)); [ $MISS -le 5 ] && echo "no reads: $S"; }
done
echo "stems with nothing in either dir: $MISS / 89"

echo ""
echo "=== 6. mapper modules available ==="
module -t avail 2>&1 | grep -iE 'bwa|samtools|minimap' | head -12

echo ""
echo "=== 7. pfam, timestamped ==="
date
ls $BASE/results/pfam_hs | grep -c '\.done$'
