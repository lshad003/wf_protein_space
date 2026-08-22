#!/bin/bash
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
WPS=/bigdata/stajichlab/lshad003/wf_protein_space

echo "== modules available =="
for M in eggnog-mapper eggnog hmmer kofamscan kofam diamond dbcan interproscan; do
  echo "--- $M ---"; module avail $M 2>&1 | grep -v "^$" | tail -3
done

echo
echo "== v1 eggnog output: how it was chunked and what format =="
ls $V1/db/eggnog_array_output | head -5
echo "total files: $(ls $V1/db/eggnog_array_output | wc -l)"
ls $V1/db/eggnog_full_output
echo "--- header of one output ---"
F=$(ls $V1/db/eggnog_array_output/* 2>/dev/null | head -1)
[ -n "$F" ] && head -6 "$F"

echo
echo "== v1 annotation scripts =="
ls $V1/pipeline/function 2>/dev/null
ls $V1/scripts | grep -i "eggnog\|kegg\|cazy\|pfam\|hmm"

echo
echo "== the split dir from v1 (95% reps were split for array jobs) =="
ls $V1/db/LsFMGC_AA_95_rep__split 2>/dev/null | head -3
echo "chunks: $(ls $V1/db/LsFMGC_AA_95_rep__split 2>/dev/null | wc -l)"

echo
echo "== candidate database locations =="
ls -d /srv/projects/db/* 2>/dev/null | head -20
ls -d /bigdata/stajichlab/shared/lib/* 2>/dev/null | head -20
ls $V1/lib 2>/dev/null | head

echo
echo "== our tiers, sizes =="
ls -lh $WPS/catalog/db/LsPS_AA_50_rep.fasta $WPS/catalog/db/LsPS_AA_95_rep.fasta
echo -n "50% reps: "; grep -c "^>" $WPS/catalog/db/LsPS_AA_50_rep.fasta
