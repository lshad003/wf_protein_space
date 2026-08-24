#!/bin/bash
#SBATCH -p stajichlab -c 8 --mem 96G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16j.log
# Step 16j: length table, then bwa-mem2 index of the supported-CDS reference.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
REF=$WPS/catalog/db/LsPS_CDS_95_supported.fasta
SUM=$WPS/results/supported_cds_extract.summary.txt

MISS=$(awk -F'\t' '$1=="missing"{print $2}' $SUM 2>/dev/null)
if [ ! -s "$REF" ] || [ "$MISS" != "0" ]; then
  echo "GUARD FAIL: ref exists=$([ -s "$REF" ] && echo yes || echo no) missing=$MISS need 0"
  exit 1
fi

/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  $WPS/scripts/step16j_lengths.py || exit 1

module load bwa-mem2 2>/dev/null
if ! command -v bwa-mem2 > /dev/null 2>&1; then
  echo "bwa-mem2 not on PATH after module load, candidates:"
  module -t avail 2>&1 | grep -i bwa
  exit 1
fi
echo "mapper: $(command -v bwa-mem2)"
bwa-mem2 version 2>&1 | head -2

bwa-mem2 index $REF
echo "=== index files ==="
ls -lh $WPS/catalog/db | grep LsPS_CDS_95_supported
echo "DONE step16j"
