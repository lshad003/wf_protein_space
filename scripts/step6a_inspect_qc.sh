#!/bin/bash
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
S=UHM574.41000

echo "== stats dirs for $S =="
ls -l $FEC/results/$S/raw_reads_stats $FEC/results/$S/clean_reads_stats
echo
echo "== contents of first file in each =="
for D in raw_reads_stats clean_reads_stats; do
  F=$(ls $FEC/results/$S/$D | head -1)
  echo "--- $D/$F ---"
  head -5 "$FEC/results/$S/$D/$F"
done
echo
echo "== qc dir =="
ls $FEC/results/$S/qc | head
echo
echo "== batch encoded in samples.csv (run directory) =="
head -1 $FEC/samples.csv
grep -E "UHM20.35743|UHM554.40991|UHM574.41000|UHM585.41009" $FEC/samples.csv
echo
echo "== distinct run directories across all samples =="
tail -n +2 $FEC/samples.csv | cut -d, -f2 | cut -d/ -f1 | sort | uniq -c | sort -rn
echo
echo "== existing WF22 QC tables in v1 =="
head -2 $V1/metagenome_read_qc_summary.tsv 2>/dev/null
echo "---"
head -2 $V1/Supplemental_Combined_Read_QC_Mapping_Metadata_updated.tsv 2>/dev/null
