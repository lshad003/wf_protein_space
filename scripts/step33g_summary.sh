#!/bin/bash
#SBATCH -J step33g_summary
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem 48G
#SBATCH -t 6:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33g_summary.%j.log
# Step 33g summary: scripts/step33g_summary.py; exits with its exit code.
W=/bigdata/stajichlab/lshad003/wf_protein_space
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step33g_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
