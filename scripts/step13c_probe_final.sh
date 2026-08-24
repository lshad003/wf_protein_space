#!/bin/bash
# Step 13c: final db probes before the GU/EU search; DUF-only x eggNOG split.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
PDB=/srv/projects/db

echo "=== 1. DUF-only supported K reps, by eggNOG status ==="
awk -F'\t' 'NR>1 && $2==1 && $6=="K" && $5==0 {if($3==1)e++; else p++}
  END{printf "dufonly with eggnog %d, dufonly pfam-only %d, total %d\n", e, p, e+p}' \
  $BASE/results/rep_classes.tsv

echo ""
echo "=== 2. Uniprot/2025_03 contents ==="
ls -lh $PDB/Uniprot/2025_03 2>/dev/null

echo ""
echo "=== 3. ncbi subdirs ==="
for D in refseq diamond mmseqs; do
  echo "--- $PDB/ncbi/$D ---"; ls -lh $PDB/ncbi/$D 2>/dev/null | head -8
done

echo ""
echo "=== 4. kaiju_nr interiors, hunting source faa ==="
ls -lh /bigdata/stajichlab/shared/db/kaiju/20260128/kaiju_nr 2>/dev/null | head -8
ls -lh /bigdata/stajichlab/shared/db/kaiju/20260128/kaiju_nr_cluster 2>/dev/null | head -8

echo ""
echo "=== 5. diamond module ==="
module -t avail 2>&1 | grep -i diamond | head -5

echo ""
echo "=== 6. array progress ==="
sacct -j 27716160 -n -X --format=State | sort | uniq -c
cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | awk 'NF==4' | sort -u | sort -t$'\t' -k4,4g
