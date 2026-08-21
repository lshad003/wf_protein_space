#!/bin/bash
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
META=/bigdata/stajichlab/lshad003/wf_protein_space/metadata
DEP=/bigdata/stajichlab/lshad003/ncbi-deposit/results
OUT=/bigdata/stajichlab/lshad003/wf_protein_space/results/step0_by_assembly_check.tsv

mkdir -p "$META"
for y in wf22 wf23 wf24; do
  cp "$DEP/${y}_stems.txt" "$META/${y}_stems.txt"
done
wc -l "$META"/wf22_stems.txt "$META"/wf23_stems.txt "$META"/wf24_stems.txt

echo -e "stem\tyear\taa_file\tcds_file" > "$OUT"
for y in wf23 wf24; do
  YR=$(echo $y | tr a-z A-Z)
  while read s; do
    aa="MISSING"; cds="MISSING"
    [ -s "$V1/by_assembly/$s.aa.fa.gz" ] && aa="present"
    [ -s "$V1/by_assembly/$s.cds.fa.gz" ] && cds="present"
    echo -e "$s\t$YR\t$aa\t$cds" >> "$OUT"
  done < "$META/${y}_stems.txt"
done
echo "== results =="
echo "rows (expect 46 with header): $(wc -l < "$OUT")"
echo "MISSING count: $(tail -n +2 "$OUT" | grep -c MISSING || true)"
echo "== any missing stems listed =="
tail -n +2 "$OUT" | awk -F'\t' '$3=="MISSING" || $4=="MISSING"'
echo "== sanity: header format of one new-stem aa file if present =="
F="$V1/by_assembly/UHM574.41000.aa.fa.gz"
[ -s "$F" ] && zcat "$F" | head -2
