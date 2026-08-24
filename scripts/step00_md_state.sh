#!/bin/bash
# After the mangled paste: what actually landed in CHATINDEX.md / PROJECT_LOG.md.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
CI=$BASE/CHATINDEX.md
PL=$BASE/PROJECT_LOG.md

echo "=== 1. unique-string counts in CHATINDEX, want 1 1 1 0 0 ==="
echo "header:      $(grep -c 'union verified, rep classes' $CI)"
echo "glob line:   $(grep -c 'glob-selected' $CI)"
echo "rule line:   $(grep -c 'boundary-defining databases' $CI)"
echo "stray echo:  $(grep -c 'PROJECT_LOG.md' $CI)"
echo "literal EOF: $(grep -cx 'EOF' $CI)"

echo ""
echo "=== 2. CHATINDEX from the afternoon header to end of file ==="
awk '/^## 2026-08-23 \(afternoon\)/{f=1} f' $CI

echo ""
echo "=== 3. PROJECT_LOG tail and duplicate check, want 1 ==="
tail -4 $PL
echo "count of today's line: $(grep -c 'union + rep classes' $PL)"

echo ""
echo "=== 4. jobs: 13f staging and progress, array ==="
tail -2 $BASE/logs/step13f.log 2>/dev/null
tail -2 $BASE/logs/step13f.diamond.err 2>/dev/null
sacct -j 27716160 -n -X --format=State | sort | uniq -c
echo "finished summaries: $(cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | awk 'NF==4' | sort -u | wc -l)"
