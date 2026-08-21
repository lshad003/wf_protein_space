#!/bin/bash
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
mkdir -p $WPS/catalog/by_assembly_new $WPS/logs $WPS/tmp

cat > $WPS/metadata/missing14_samples.csv << 'CSV'
stem,shotgun
UHM585.41009,NA
UHM617.41023,NA
UHM624.41024,NA
UHM626.41025,NA
UHM628.41026,NA
UHM632.41027,NA
UHM633.41028,NA
UHM634.41029,NA
UHM636.41030,NA
UHM644.41031,NA
UHM648.41909,NA
UHM649.41032,NA
UHM656.41033,NA
UHM666.41034,NA
CSV

echo "csv data rows (expect 14): $(tail -n +2 $WPS/metadata/missing14_samples.csv | wc -l)"
echo
echo "=== contig file check ==="
MISS=0
tail -n +2 $WPS/metadata/missing14_samples.csv | while IFS=, read S REST; do
  F=$FEC/results/$S/scaffolds/${S}_R.fa.gz
  if [ -s "$F" ]; then
    echo "OK   $S  $(stat -c %s "$F") bytes"
  else
    echo "MISS $S  $F"
  fi
done
echo
echo "=== prodigal module ==="
module load prodigal && prodigal -v 2>&1 | head -2
