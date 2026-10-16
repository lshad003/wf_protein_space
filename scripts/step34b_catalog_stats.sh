#!/bin/bash
#SBATCH -J step34b_catalog_stats
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 6:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34b_catalog_stats.%j.log
# Step 34b: scripts/step34b_catalog_stats.py (see its header); exits with its exit code.
W=/bigdata/stajichlab/lshad003/wf_protein_space
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step34b_catalog_stats.py
rc=$?
echo "exit $rc"; date
exit $rc
