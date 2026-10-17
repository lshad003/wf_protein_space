#!/bin/bash
#SBATCH -J step34c3b
#SBATCH -p stajichlab -c 2 --mem 48G --time=6:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34c3b_checks.%j.log
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34c3b_checks.py
echo "exit $?"; date
