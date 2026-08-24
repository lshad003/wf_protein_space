#!/bin/bash
# Step 13g: watch the nr_cluster screen and the mapping array.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
echo "=== queue ==="
squeue -u lshad003 -o "%.10i %.9P %.14j %.2t %.10M %R" | head -14
echo ""
echo "=== 13f log and diamond progress ==="
tail -3 $BASE/logs/step13f.log 2>/dev/null
tail -2 $BASE/logs/step13f.diamond.err 2>/dev/null
echo ""
echo "=== array fractions, low to high ==="
sacct -j 27716160 -n -X --format=State | sort | uniq -c
X=$(cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | awk 'NF==4' | sort -u | sort -t$'\t' -k4,4g)
echo "$X" | head -5; echo "..."; echo "$X" | tail -2; echo "finished: $(echo "$X" | wc -l)"
