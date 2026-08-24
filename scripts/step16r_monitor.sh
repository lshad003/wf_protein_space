#!/bin/bash
# Step 16r: afternoon monitor. Union numbers, family rollup, array progress.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space

echo "=== 1. job states ==="
sacct -j 27716158,27716159 -n -X --format=JobID%12,JobName%16,State%12,ExitCode,Elapsed
sacct -j 27716160 -n -X --format=State | sort | uniq -c

echo ""
echo "=== 2. step12d union table ==="
tail -4 $BASE/logs/step12d.log 2>/dev/null
[ -s $BASE/results/union_rep_level.tsv ] && cat $BASE/results/union_rep_level.tsv

echo ""
echo "=== 3. step12e family union ==="
tail -6 $BASE/logs/step12e.log 2>/dev/null

echo ""
echo "=== 4. array progress, lowest and highest mapped fractions ==="
N=$(ls $BASE/results/mapping_summary 2>/dev/null | wc -l)
echo "summaries: $N / 89"
cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | sort -t$'\t' -k4,4g | head -8
echo "..."
cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | sort -t$'\t' -k4,4g | tail -3

echo ""
echo "=== 5. problems in array logs, if any ==="
grep -l -E 'GUARD FAIL|MISSING READS|no row' $BASE/logs/step16q.*.log 2>/dev/null | head -5
