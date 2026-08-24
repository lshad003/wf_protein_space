#!/bin/bash
# Step 16d: inventory of renamed inputs, WF22 CDS sources, reads, emapper tmp.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
OLD=$FEC/Proteins/by_assembly

echo "=== 1. step16b state and log tail ==="
sacct -j 27713254 -n -X --format=State,Elapsed,MaxRSS
tail -8 $BASE/logs/step16b.log 2>/dev/null

echo ""
echo "=== 2. where step2 took its catalog input from ==="
grep -nE 'input_pep|by_assembly|INPUT|cat .*fasta|\.fa' $BASE/scripts/step2_cluster.sh | head -10

echo ""
echo "=== 3. renamed files per cohort in input_pep and input_cds ==="
for D in input_pep input_cds; do
  W22=0; W23=0; W24=0
  for S in $(cat $BASE/metadata/wf22_stems_44.txt); do [ -s $BASE/catalog/$D/$S.fasta ] && W22=$((W22+1)); done
  for S in $(cat $BASE/metadata/wf23_stems.txt); do [ -s $BASE/catalog/$D/$S.fasta ] && W23=$((W23+1)); done
  for S in $(cat $BASE/metadata/wf24_stems.txt); do [ -s $BASE/catalog/$D/$S.fasta ] && W24=$((W24+1)); done
  echo "$D: total $(ls $BASE/catalog/$D 2>/dev/null | wc -l), wf22 $W22/44, wf23 $W23/9, wf24 $W24/36"
done

echo ""
echo "=== 4. WF22 cds sources in shared by_assembly ==="
M=""
for S in $(cat $BASE/metadata/wf22_stems_44.txt); do
  [ -s $OLD/$S.cds.fa.gz ] || M="$M $S"
done
echo "wf22 stems lacking cds.fa.gz: $(echo $M | wc -w)/44"
echo $M | tr ' ' '\n' | head -5

echo ""
echo "=== 5. spot check, one stem: pep vs cds sequence counts ==="
S=$(head -1 $BASE/metadata/wf23_stems.txt)
echo "stem: $S"
grep -c '>' $BASE/catalog/input_pep/$S.fasta
grep -c '>' $BASE/catalog/input_cds/$S.fasta

echo ""
echo "=== 6. reads: inside one qc dir ==="
S=$(head -1 $BASE/metadata/wf22_stems_44.txt)
ls -l $FEC/results/$S/qc 2>/dev/null | head -12

echo ""
echo "=== 7. emapper tmp: count, size, and eggnog2 completeness ==="
ls -d $BASE/emappertmp_dmdn_* 2>/dev/null | wc -l
du -shc $BASE/emappertmp_dmdn_* 2>/dev/null | tail -1
ls $BASE/results/eggnog2/*.emapper.annotations 2>/dev/null | wc -l

echo ""
echo "=== 8. Pfam progress, timestamped ==="
date
ls $BASE/results/pfam_hs | grep -c '\.done$'
