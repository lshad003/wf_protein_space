#!/bin/bash
# Step 16o: morning verdicts. Pfam completion, 12d union, 16m pilot. Read, not guess.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space

echo "=== 1. final states: pfam, union 27713842, pilot 27713872 ==="
sacct -j 27712560,27712562,27713842,27713872 -S 2026-08-22 -n -X \
  --format=JobID%16,JobName%14,State%12,ExitCode,Elapsed,End | head -20

echo ""
echo "=== 2. pfam done count, want 619 ==="
ls $BASE/results/pfam_hs | grep -c '\.done$'

echo ""
echo "=== 3. step12d union log, want table + cross-check pass ==="
tail -15 $BASE/logs/step12d.log 2>/dev/null

echo ""
echo "=== 4. union outputs on disk ==="
ls -l $BASE/results/union_rep_level.tsv $BASE/results/known_union_ids.txt \
      $BASE/results/pfam_hit_ids.txt $BASE/results/eggnog_hit_ids.txt 2>/dev/null

echo ""
echo "=== 5. pilot log, want summary line + nonzero gene count ==="
tail -12 $BASE/logs/step16m_pilot.log 2>/dev/null

echo ""
echo "=== 6. pilot outputs ==="
ls -l $BASE/results/mapping_summary $BASE/results/counts 2>/dev/null

echo ""
echo "=== 7. pilot bwa stderr tail, mapping health ==="
tail -5 $BASE/logs/step16m_pilot.bwa.err 2>/dev/null
