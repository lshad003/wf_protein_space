#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 8G --time 6:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16h.log
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
N=$(ls $WPS/catalog/input_cds | wc -l)
OK=$(grep -c 'ID SET IDENTICAL' $WPS/logs/step16e.log 2>/dev/null)
if [ "$N" != "89" ] || [ "$OK" != "1" ]; then
  echo "GUARD FAIL: input_cds=$N need 89, step16e verified=$OK need 1"
  exit 1
fi
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  $WPS/scripts/step16h_supported_cds.py
echo "headers in reference: $(grep -c '>' $WPS/catalog/db/LsPS_CDS_95_supported.fasta), expect 2069453"
