#!/bin/bash
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
NEW=$WPS/catalog/by_assembly_new
OUT=$WPS/results/step1_protein_counts_45.tsv

echo -e "stem\twf_year\tprediction_batch\tn_proteins\tcontig_bytes" > $OUT
while read S; do
  if grep -qx "$S" $WPS/metadata/wf23_stems.txt; then YR=WF23; else YR=WF24; fi
  if [ -s "$NEW/$S.aa.fa.gz" ]; then B=new_2026; else B=preexisting; fi
  N=$(grep -c "^>" $WPS/catalog/input_pep/$S.fasta)
  C=$(stat -c %s "$FEC/results/$S/scaffolds/${S}_R.fa.gz" 2>/dev/null || echo NA)
  echo -e "$S\t$YR\t$B\t$N\t$C" >> $OUT
done < $WPS/metadata/new45_stems.txt

echo "rows (expect 46): $(wc -l < $OUT)"
echo "== totals by batch =="
awk -F'\t' 'NR>1 {n[$3]++; s[$3]+=$4} END {for (b in n) printf "%s: %d samples, %d proteins, mean %d\n", b, n[b], s[b], s[b]/n[b]}' $OUT
echo "== totals by year =="
awk -F'\t' 'NR>1 {n[$2]++; s[$2]+=$4} END {for (y in n) printf "%s: %d samples, %d proteins, mean %d\n", y, n[y], s[y], s[y]/n[y]}' $OUT
