#!/bin/bash
#SBATCH -J step34a_origin
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem 64G
#SBATCH -t 6:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34a_origin.%j.log
# Step 34a origin: scripts/step34a_origin.py (see its header); exits with its exit code.
W=/bigdata/stajichlab/lshad003/wf_protein_space
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step34a_origin.py
rc=$?
echo "exit $rc"; date
exit $rc
