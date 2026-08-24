#!/bin/bash
# Step 20a: job checks + free GU taxonomy from hit titles + existing scaffold tax.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
export LC_ALL=C

echo "=== 1. job states ==="
sacct -j 27721791,27721810 -n -X --format=JobID%12,JobName%16,State%12,ExitCode,Elapsed
echo "--- 17e ---"; tail -6 $WPS/logs/step17e.log 2>/dev/null
echo "--- 13i ---"; tail -12 $WPS/logs/step13i.log 2>/dev/null

echo ""
echo "=== 2. GU taxonomy free from nr hit titles (top organisms) ==="
awk -F'\t' '$6<=1e-10 && $8>=50 && $9>=50 {print $1"\t"$10}' \
  $WPS/results/step13_nrclust20260128_hits.tsv \
| awk -F'\t' '{if (match($2,/\[[^]]+\]/)) {s=substr($2,RSTART+1,RLENGTH-2); if(!seen[$1]++) print s}}' \
| sort | uniq -c | sort -rn | head -15

echo ""
echo "=== 3. existing scaffold taxonomy in Fecal ==="
ls $FEC/results_scaffold_classify_mmseqs 2>/dev/null | head -5
ls $FEC/results_kraken2 2>/dev/null | head -5
S=$(head -1 $WPS/metadata/wf22_stems_44.txt)
echo "--- files for $S ---"
ls -l $FEC/results_scaffold_classify_mmseqs/${S}* 2>/dev/null | head -4
ls -l $FEC/results_kraken2/${S}* 2>/dev/null | head -4
