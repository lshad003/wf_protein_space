#!/bin/bash
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
META=/bigdata/stajichlab/lshad003/wf_protein_space/metadata
OUT=/bigdata/stajichlab/lshad003/wf_protein_space/results/step0_by_assembly_check.tsv

echo "== by_assembly dir overview =="
ls "$V1/by_assembly" | head -20
echo "...total files:"
ls "$V1/by_assembly" | wc -l
echo
echo "== check 45 new stems for existing aa/cds predictions =="
echo -e "stem\taa_file\tcds_file" > "$OUT"
for f in "$META/wf23_stems.txt" "$META/wf24_stems.txt"; do
  while read s; do
    aa="MISSING"; cds="MISSING"
    [ -s "$V1/by_assembly/$s.aa.fa.gz" ] && aa="$V1/by_assembly/$s.aa.fa.gz"
    [ -s "$V1/by_assembly/$s.cds.fa.gz" ] && cds="$V1/by_assembly/$s.cds.fa.gz"
    echo -e "$s\t$aa\t$cds" >> "$OUT"
  done < "$f"
done
echo "wrote $OUT"
echo "MISSING aa count: $(grep -c MISSING <(cut -f2 "$OUT" | tail -n +2))"
echo
echo "== pipeline scripts to adopt =="
ls "$V1/pipeline"
echo
echo "== gene_prediction dir (alternate location) =="
ls "$V1/gene_prediction" | head -20
