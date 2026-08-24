#!/bin/bash
# Step 16s: array progress (non-empty only), stranded partials, 12f verdict,
# genome-database recon for the GU/EU search.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal

echo "=== 1. array task states ==="
sacct -j 27716160 -n -X --format=State | sort | uniq -c

echo ""
echo "=== 2. finished summaries (non-empty), spread of mapped fraction ==="
ls -l $BASE/results/mapping_summary | awk '$5>0 && NF>2 {print $NF}' | wc -l
cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | awk 'NF==4' \
  | sort -t$'\t' -k4,4g | awk 'NR<=5{print} END{print "...highest:"; print}'

echo ""
echo "=== 3. stranded uncompressed counts, want none ==="
ls $BASE/results/counts 2>/dev/null | grep -v '\.gz$' | head -5

echo ""
echo "=== 4. step12f verdict ==="
tail -5 $BASE/logs/step12f.log 2>/dev/null

echo ""
echo "=== 5. genome DB candidates on cluster ==="
echo "--- Fecal/gtdbtk_db ---"
ls $FEC/gtdbtk_db 2>/dev/null | head -8
echo "--- shared db dirs ---"
ls /bigdata/stajichlab/shared/db 2>/dev/null | head -10
ls /srv/projects/db 2>/dev/null | head -10
echo "--- lshad003 databases ---"
ls /bigdata/stajichlab/lshad003/ncbi-deposit/wf_databases 2>/dev/null | head -6
