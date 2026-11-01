#!/bin/bash
#SBATCH -J step36b_nearmiss
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 2:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step36b_nearmiss.%j.log
# Step 36b B: near-miss table from the hmmsearch --max array (scripts/step36b_nearmiss.py).
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step36b_nearmiss.py
rc=$?
echo "exit $rc"; date
exit $rc
