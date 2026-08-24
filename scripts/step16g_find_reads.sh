#!/bin/bash
# Step 16g: step16e verdict, then locate fastqs the way v1 did.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1

echo "=== 1. step16e state and log tail ==="
sacct -j 27713782 -n -X --format=State,Elapsed
tail -12 $BASE/logs/step16e.log 2>/dev/null

echo ""
echo "=== 2. Fecal top level ==="
ls $FEC

echo ""
echo "=== 3. probe common read dirs for one stem ==="
S=UHM102.10840
for D in input reads data fastq filtered clean; do
  ls $FEC/$D/${S}* 2>/dev/null | head -3
done

echo ""
echo "=== 4. v1 layout and its mapping conventions ==="
ls $V1
echo "--- v1 scripts mentioning fastq or mapping ---"
grep -lE 'fq\.gz|fastq|bwa|bowtie|minimap' $V1/scripts/*.sh 2>/dev/null | head -6
for F in $(grep -lE 'fq\.gz|fastq|bwa|bowtie|minimap' $V1/scripts/*.sh 2>/dev/null | head -2); do
  echo "--- $F ---"
  grep -nE 'fq\.gz|fastq|bwa|bowtie|minimap|module load|R1|R2' $F | head -10
done

echo ""
echo "=== 5. pfam progress, timestamped ==="
date
ls $BASE/results/pfam_hs | grep -c '\.done$'
