#!/bin/bash
#SBATCH -J step34c2_genus
#SBATCH -p stajichlab -c 2 --mem 48G --time=12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34c2_genus.%j.log
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34c2_genus.py
echo "exit $?"; date
