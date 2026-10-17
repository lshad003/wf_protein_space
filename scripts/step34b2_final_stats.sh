#!/bin/bash
#SBATCH -J step34b2
#SBATCH -p stajichlab -c 2 --mem 48G --time=6:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34b2_final_stats.%j.log
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34b2_final_stats.py
echo "exit $?"; date
