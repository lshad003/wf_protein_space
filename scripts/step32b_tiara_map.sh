#!/bin/bash
#SBATCH -J step32b_tiara_map
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 4:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step32b_tiara_map.%j.log
# Step 32b: scripts/step32b_tiara_map.py (see its header); exits with its exit code.
W=/bigdata/stajichlab/lshad003/wf_protein_space
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step32b_tiara_map.py
rc=$?
echo "exit $rc"; date
exit $rc
