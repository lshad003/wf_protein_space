#!/bin/bash
# Step 17a: one sample table for all 89, with cohort and treatment.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
OUT=$WPS/metadata/samples_89_treatment.tsv
export LC_ALL=C

printf "stem\tcohort\ttreatment\tegg_mass\n" > $OUT

# WF22 from v1 metadata, subset to our 44
awk -F'\t' 'NR>1 && $3!="" {print $1"\t"$3"\t"$5}' $V1/last_wood_frog_2022_sample_metadata.tsv | sort > /tmp/w22.tsv
sort $WPS/metadata/wf22_stems_44.txt > /tmp/w22ids.txt
join -t$'\t' /tmp/w22ids.txt /tmp/w22.tsv | awk -F'\t' '{print $1"\tWF22\t"$2"\t"$3}' >> $OUT

# WF23 and WF24 from project metadata
awk -F'\t' 'NR>1 {print $1"\t"$3"\t"$5"\t"$4}' $WPS/metadata/wf23_wf24_sequenced_metadata.tsv >> $OUT

echo "rows: $(($(wc -l < $OUT) - 1)), expect 89"
echo ""
echo "=== treatment counts per cohort ==="
awk -F'\t' 'NR>1{print $2"\t"$3}' $OUT | sort | uniq -c
echo ""
echo "=== missing from table ==="
cut -f1 $OUT | tail -n +2 | sort > /tmp/have.txt
sort $WPS/results/cohort_map.tsv | cut -f1 | sort > /tmp/want.txt
comm -13 /tmp/have.txt /tmp/want.txt
