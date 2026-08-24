#!/bin/bash
# Step 16k: join 89 stems to read_manifest.csv, verify every R1/R2 exists.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
MAN=$FEC/read_manifest.csv
OUT=$BASE/results/reads_89.tsv

echo "=== 1. step16j state and log tail ==="
sacct -j 27713810 -n -X --format=State,Elapsed,MaxRSS
tail -8 $BASE/logs/step16j.log 2>/dev/null

echo ""
echo "=== 2. manifest join, all 89 ==="
NOMAN=0; NOFILE=0; N=0
printf "stem\tR1\tR2\n" > $OUT
for S in $(cat $BASE/metadata/wf22_stems_44.txt $BASE/metadata/wf23_stems.txt $BASE/metadata/wf24_stems.txt); do
  N=$((N+1))
  LINE=$(awk -F, -v s="$S" '$1==s{print; exit}' $MAN)
  if [ -z "$LINE" ]; then NOMAN=$((NOMAN+1)); [ $NOMAN -le 10 ] && echo "not in manifest: $S"; continue; fi
  R1=$(echo "$LINE" | cut -d, -f2)
  R2=$(echo "$LINE" | cut -d, -f3)
  if [ -s "$R1" ] && [ -s "$R2" ]; then
    printf "%s\t%s\t%s\n" "$S" "$R1" "$R2" >> $OUT
  else
    NOFILE=$((NOFILE+1)); [ $NOFILE -le 10 ] && echo "manifest ok, file missing: $S"
  fi
done
echo "stems checked: $N, not in manifest: $NOMAN, files missing: $NOFILE"
echo "reads_89.tsv rows: $(($(wc -l < $OUT) - 1)), expect 89"

echo ""
echo "=== 3. one manifest line per cohort, filename patterns ==="
for S in $(head -1 $BASE/metadata/wf22_stems_44.txt) $(head -1 $BASE/metadata/wf23_stems.txt) $(head -1 $BASE/metadata/wf24_stems.txt); do
  awk -F, -v s="$S" '$1==s{print; exit}' $MAN
done

echo ""
echo "=== 4. step11 family-join logic, for the union rollup ==="
grep -nE 'clusters_50|clusters_95|family|rep' $BASE/scripts/step11_interaction_and_supported.sh | head -15

echo ""
echo "=== 5. pfam, timestamped ==="
date
ls $BASE/results/pfam_hs | grep -c '\.done$'
