#!/bin/bash
#SBATCH -p stajichlab -c 2 --mem 4G --time 2:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16e.log
# Step 16e: rename 44 WF22 CDS files with the exact step1d transform, then
# verify ID sets against v1 input_pep, the literal source step2 catted.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
OLD=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly
V1PEP=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1/input_pep
CDS=$WPS/catalog/input_cds

echo "=== header formats, first stem ==="
S=$(head -1 $WPS/metadata/wf22_stems_44.txt)
echo "--- v1 pep ---";  grep '>' $V1PEP/$S.fasta | head -2
echo "--- raw cds ---"; pigz -dc $OLD/$S.cds.fa.gz | grep '>' | head -2

BAD=0
while read S; do
  export BIOSAMPLE=$S
  pigz -dc $OLD/$S.cds.fa.gz | perl -p -e 's/^>(\S+).+/>$ENV{BIOSAMPLE}__$1/' > $CDS/$S.fasta
  NP=$(grep -c '>' $V1PEP/$S.fasta)
  NC=$(grep -c '>' $CDS/$S.fasta)
  [ "$NP" = "$NC" ] || { echo "COUNT MISMATCH $S pep=$NP cds=$NC"; BAD=$((BAD+1)); }
done < $WPS/metadata/wf22_stems_44.txt
echo "=== stems with count mismatch: $BAD of 44, expect 0 ==="

echo "=== exact ID set, first stem, first-word compare ==="
S=$(head -1 $WPS/metadata/wf22_stems_44.txt)
A=$(grep '>' $V1PEP/$S.fasta | awk '{print $1}' | sort | md5sum | awk '{print $1}')
B=$(grep '>' $CDS/$S.fasta   | awk '{print $1}' | sort | md5sum | awk '{print $1}')
echo "pep $A"; echo "cds $B"
[ "$A" = "$B" ] && echo "ID SET IDENTICAL" || echo "ID SET DIFFERS, STOP"
echo "=== input_cds files now: $(ls $CDS | wc -l), expect 89 ==="
