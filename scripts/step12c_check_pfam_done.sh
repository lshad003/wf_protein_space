#!/bin/bash
# Step 12c: verify Pfam hmmsearch completion before any Pfam number is used,
# and locate the inputs needed for the eggNOG+Pfam union.

BASE=/bigdata/stajichlab/lshad003/wf_protein_space
PF=$BASE/results/pfam_hs

echo "=== 1. job state histogram, jobs 27712560 27712562 ==="
sacct -j 27712560,27712562 -S 2026-08-20 -n -X --format=State | sort | uniq -c
echo "--- non-completed tasks, if any (top 20) ---"
sacct -j 27712560,27712562 -S 2026-08-20 -n -X --format=JobID%22,State%12,Elapsed | grep -vE 'COMPLETED|RUNNING|PENDING' | head -20

echo ""
echo "=== 2. .done files, expect 619 ==="
ls $PF 2>/dev/null | grep -c '\.done$'

echo ""
echo "=== 3. pfam_hs contents by extension ==="
ls $PF 2>/dev/null | awk -F. 'NF>1 {print $NF}' | sort | uniq -c | sort -rn

echo ""
echo "=== 4. zero-size non-done files, expect 0 ==="
ls -l $PF 2>/dev/null | awk '$5==0 && $NF !~ /\.done$/' | wc -l
ls -l $PF 2>/dev/null | awk '$5==0 && $NF !~ /\.done$/ {print $NF}' | head -5

echo ""
echo "=== 5. output settings used by step12b ==="
grep -nE 'hmmsearch|tblout|domtbl|cut_ga|noali|-E ' $BASE/scripts/step12b_pfam_hmmsearch.sh | head -15

echo ""
echo "=== 6. first data lines of one output file ==="
SAMPLE=$(ls $PF | grep -v '\.done$' | head -1)
echo "sample: $SAMPLE"
case "$SAMPLE" in
  *.gz) zcat "$PF/$SAMPLE" | grep -v '^#' | head -3 ;;
  *)    grep -v '^#' "$PF/$SAMPLE" | head -3 ;;
esac

echo ""
echo "=== 7. eggNOG annotation source used by step10 ==="
grep -nE 'emapper|annotation' $BASE/scripts/step10_annotation_summary.sh | head -12
