#!/bin/bash
# Step 14a: locate AntiFam HMMs and hmmsearch before running the screen.
PDB=/srv/projects/db
SDB=/bigdata/stajichlab/shared/db

echo "=== 1. AntiFam anywhere in the db roots ==="
ls -d $PDB/*ntifam* $PDB/*ANTIFAM* $SDB/*ntifam* 2>/dev/null
ls $PDB/pfam 2>/dev/null | head -10
ls $PDB/PFAM 2>/dev/null | head -10

echo ""
echo "=== 2. what our Pfam step used ==="
grep -nE 'DB=|Pfam-A|hmm' /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step12b_pfam_hmmsearch.sh | head -8

echo ""
echo "=== 3. hmmer module ==="
module -t avail 2>&1 | grep -i hmmer | head -5

echo ""
echo "=== 4. query file ready? ==="
ls -lh /bigdata/stajichlab/lshad003/wf_protein_space/results/supported_U_reps.faa
grep -c '>' /bigdata/stajichlab/lshad003/wf_protein_space/results/supported_U_reps.faa
