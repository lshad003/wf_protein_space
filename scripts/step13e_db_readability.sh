#!/bin/bash
# Step 13e: which candidate databases exist and which installed diamond reads
# them; ncbi/diamond and ncbi/mmseqs contents; step13d postmortem; array.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
PDB=/srv/projects/db

echo "=== 1. step13d final state and log tail ==="
sacct -j 27716261 -n -X --format=State,ExitCode,Elapsed
tail -8 $BASE/logs/step13d.log 2>/dev/null

echo ""
echo "=== 2. diamond binaries available ==="
module load diamond 2>/dev/null; diamond version 2>&1 | head -1
ls -d /opt/linux/rocky/8.x/x86_64/pkgs/diamond/* 2>/dev/null
ls /bigdata/stajichlab/lshad003/condaenvs/*/bin/diamond 2>/dev/null

echo ""
echo "=== 3. dbinfo readability, module diamond ==="
for D in $PDB/Uniprot/uniref90.dmnd $PDB/Uniprot/2025_03/uniref100.dmnd $PDB/Uniprot/2025_03/uniprot_sprot.dmnd; do
  echo "--- $D ---"
  diamond dbinfo --db $D 2>&1 | grep -iE 'version|sequences|letters|rror' | head -4
done

echo ""
echo "=== 4. ncbi/diamond current and 20250623 ==="
ls -lh $PDB/ncbi/diamond/20260128 2>/dev/null | head -8
ls -lh $PDB/ncbi/diamond/20250623 2>/dev/null | head -8
N1=$(ls $PDB/ncbi/diamond/20260128/*.dmnd 2>/dev/null | head -1)
if [ -n "$N1" ]; then echo "--- dbinfo $N1 ---"; diamond dbinfo --db $N1 2>&1 | grep -iE 'version|sequences|letters|rror' | head -4; fi

echo ""
echo "=== 5. ncbi/mmseqs, named targets ==="
ls -lh $PDB/ncbi/mmseqs 2>/dev/null | grep -iE 'uniref|swiss|refseq|nr|uniprot' | head -12
ls $PDB/ncbi/mmseqs/archive 2>/dev/null | head -15

echo ""
echo "=== 6. kaiju group fastas, source protein sets ==="
ls -lh /bigdata/stajichlab/shared/db/kaiju/20260128/kaiju_nr/group_fastas 2>/dev/null | head -8

echo ""
echo "=== 7. array progress ==="
sacct -j 27716160 -n -X --format=State | sort | uniq -c
cat $BASE/results/mapping_summary/*.summary.tsv 2>/dev/null | awk 'NF==4' | sort -u | sort -t$'\t' -k4,4g
