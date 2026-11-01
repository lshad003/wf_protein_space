#!/bin/bash
#SBATCH -J step36b_prep
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 2:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step36b_prep.%j.log
# Step 36b prep: member list, faidx extraction, deduplication and chunking (scripts/step36b_prep.py). $1 = number of chunks.
W=/bigdata/stajichlab/lshad003/wf_protein_space; D=$W/results/step36b; PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
module load samtools/1.22.1
hostname; date
$PY $W/scripts/step36b_prep.py members || exit 1
if [ -e $D/members.faa ]; then echo "REFUSING: $D/members.faa exists"; exit 1; fi
samtools faidx $W/catalog/db/LsPS_AA_all.fasta.gz -r $D/ids.txt > $D/members.faa || { echo "exit 1 (faidx)"; exit 1; }
$PY $W/scripts/step36b_prep.py split ${1:-32}
rc=$?
echo "exit $rc"; date
exit $rc
