#!/bin/bash
#SBATCH -p epyc
#SBATCH -c 2
#SBATCH --mem 96G
#SBATCH -t 2-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_%x.%j.log
# Step 34j: run one Step 34j python script ($1, in scripts/) with rf_py39.
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/$1
rc=$?; echo "exit $rc"; date; exit $rc
