#!/bin/bash
#SBATCH -J step33c2_summary
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem 48G
#SBATCH -t 6:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33c2_summary.%j.log
# Step 33c v2 summary: scripts/step33c2_summary.py; exits with its exit code (3 = negative-control gate failed).
W=/bigdata/stajichlab/lshad003/wf_protein_space
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step33c2_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
