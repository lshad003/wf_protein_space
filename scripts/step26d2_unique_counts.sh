#!/usr/bin/bash -l
#SBATCH -p stajichlab -c 2 --mem 48gb --time=12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step26d2.%j.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step26d2_unique_counts.py
