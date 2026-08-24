#!/bin/bash
# Step 13a: full database inventory for the GU/EU search + running-job status.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
SDB=/bigdata/stajichlab/shared/db
PDB=/srv/projects/db

echo "=== 1. step12f verdict ==="
sacct -j 27716207 -n -X --format=State,ExitCode,Elapsed
tail -6 $BASE/logs/step12f.log 2>/dev/null

echo ""
echo "=== 2. array: states and finished summaries ==="
sacct -j 27716160 -n -X --format=State | sort | uniq -c
cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | awk 'NF==4' | sort -u \
  | sort -t$'\t' -k4,4g | awk '{print} END{print "finished: " NR}'

echo ""
echo "=== 3. shared/db, complete listing ==="
ls -1 $SDB

echo ""
echo "=== 4. /srv/projects/db, complete listing ==="
ls -1 $PDB

echo ""
echo "=== 5. named probes ==="
for D in $SDB/genomes $SDB/kaiju $PDB/alphafold; do
  echo "--- $D ---"
  ls $D 2>/dev/null | head -12
done
echo "--- alphafold/uniref90 ---"
ls -lh $PDB/alphafold/uniref90 2>/dev/null | head -4
echo "--- anything gtdb/uniref/refseq/nr named, both roots ---"
ls -d $SDB/*gtdb* $SDB/*GTDB* $SDB/*niref* $SDB/*efseq* $SDB/nr* \
      $PDB/*gtdb* $PDB/*GTDB* $PDB/*niref* $PDB/*efseq* $PDB/nr* 2>/dev/null
