#!/bin/bash
# Step 16l: 12d v1 verdict, 16j progress, read-path uniqueness, ref size, pfam.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space

echo "=== 1. step12d v1, job 27713815, expect GUARD FAIL ==="
sacct -j 27713815 -n -X --format=State,ExitCode,Elapsed
tail -3 $BASE/logs/step12d.log 2>/dev/null

echo ""
echo "=== 2. step16j, job 27713810 ==="
sacct -j 27713810 -n -X --format=State,Elapsed,MaxRSS
tail -4 $BASE/logs/step16j.log 2>/dev/null

echo ""
echo "=== 3. read path uniqueness across 89, expect 89 and 89, no dups ==="
echo "distinct R1: $(cut -f2 $BASE/results/reads_89.tsv | tail -n +2 | sort -u | wc -l)"
echo "distinct R2: $(cut -f3 $BASE/results/reads_89.tsv | tail -n +2 | sort -u | wc -l)"
cut -f2 $BASE/results/reads_89.tsv | tail -n +2 | sort | uniq -d | head -5

echo ""
echo "=== 4. supported CDS reference size, from the length table ==="
awk -F'\t' '{s+=$2} END{printf "seqs %d, total bp %d, mean %.1f\n", NR, s, s/NR}' \
  $BASE/results/supported_cds_lengths.tsv 2>/dev/null

echo ""
echo "=== 5. pfam, timestamped ==="
date
ls $BASE/results/pfam_hs | grep -c '\.done$'
