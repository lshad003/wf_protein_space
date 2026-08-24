#!/bin/bash
# Step 16a: preflight for abundance mapping while Pfam runs.
# Verifies the facts the mapping design depends on. No heavy work.

BASE=/bigdata/stajichlab/lshad003/wf_protein_space

echo "=== 1. what is running now, partitions Pfam occupies ==="
squeue -u lshad003 -o "%.10i %.9P %.12j %.2t %.10M %.5C %R" | head -15

echo ""
echo "=== 2. supported cluster set: where step4 wrote it ==="
grep -nE 'results/|\.tsv|\.txt' $BASE/scripts/step4_nonsingleton.sh | head -12
ls -l $BASE/results | grep -iE 'support|singleton|step4|step3'

echo ""
echo "=== 3. cluster membership table: where step3 read it ==="
grep -nE 'createtsv|cluster.*tsv|catalog/' $BASE/scripts/step3_cluster_composition.sh | head -10

echo ""
echo "=== 4. nucleotide gene sequences: do they exist and where ==="
echo "--- step1c prodigal outputs, the 14 predicted here ---"
grep -nE '\-d |\.fna|\.ffn|prodigal' $BASE/scripts/step1c_prodigal_array.sh | head -10
echo "--- step1d: did renaming cover nucleotide files ---"
grep -nE '\.fna|\.ffn|\.faa' $BASE/scripts/step1d_rename_45.sh | head -10
echo "--- step0b: where pre-existing predictions live ---"
grep -nE '/bigdata|\.fna|\.ffn|\.faa' $BASE/scripts/step0b_check_by_assembly.sh | head -10

echo ""
echo "=== 5. reads for the 89: where step6a found them ==="
grep -nE '/bigdata|fastq|fq\.gz|reads' $BASE/scripts/step6a_inspect_qc.sh | head -10

echo ""
echo "=== 6. metadata schema, for the WF24 treatment count check ==="
head -1 $BASE/metadata/wf23_wf24_sequenced_metadata.tsv

echo ""
echo "=== 7. disk headroom on /bigdata ==="
df -h /bigdata
