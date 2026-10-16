#!/bin/bash
#SBATCH -J step33e_tiara_vs_kraken
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 4:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33e_tiara_vs_kraken.%j.log
# Step 33e: scripts/step33e_tiara_vs_kraken.py (see its header); exits with its exit code.
W=/bigdata/stajichlab/lshad003/wf_protein_space
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step33e_tiara_vs_kraken.py
rc=$?
echo "exit $rc"; date
exit $rc
