#!/bin/bash
# Step 13b: strict-K at supported level; targeted probes for uniref90, ncbi,
# kaiju source fasta, uniprot; array progress.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
SDB=/bigdata/stajichlab/shared/db
PDB=/srv/projects/db

echo "=== 1. supported-level strict-K from rep_classes.tsv ==="
awk -F'\t' 'NR>1 && $2==1 {k[$6]++; if($6=="K" && $5==0) d++}
  END{printf "supported K %d KWP %d U %d | K dufonly %d | strict K %d, strict GUside %d\n",
      k["K"],k["KWP"],k["U"],d,k["K"]-d,k["U"]+d}' $BASE/results/rep_classes.tsv

echo ""
echo "=== 2. uniref90 inside alphafold versions ==="
for V in 2.1.2 2.3.0 3.0.0 3.0.0_11122024; do
  ls -lh $PDB/alphafold/$V/uniref90/uniref90.fasta 2>/dev/null
done
echo "--- alphafold/2.3.0 contents ---"
ls $PDB/alphafold/2.3.0 2>/dev/null | head -12

echo ""
echo "=== 3. ncbi dirs, both roots ==="
for D in $SDB/ncbi $PDB/ncbi $PDB/NCBI; do
  echo "--- $D ---"; ls $D 2>/dev/null | head -12
done

echo ""
echo "=== 4. kaiju builds, hunting a source .faa ==="
ls -lh $SDB/kaiju/20260128 2>/dev/null | head -10
echo "--- srv kaiju ---"
ls $PDB/kaiju 2>/dev/null | head -8

echo ""
echo "=== 5. uniprot family ==="
for D in $PDB/uniprot $PDB/Uniprot $PDB/Swissprot; do
  echo "--- $D ---"; ls -lh $D 2>/dev/null | head -8
done

echo ""
echo "=== 6. array progress, finished summaries by mapped fraction ==="
sacct -j 27716160 -n -X --format=State | sort | uniq -c
cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | awk 'NF==4' | sort -u | sort -t$'\t' -k4,4g
