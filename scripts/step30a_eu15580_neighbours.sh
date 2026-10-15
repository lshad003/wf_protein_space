#!/bin/bash
#SBATCH -J step30a_nb
#SBATCH -p epyc
#SBATCH -c 4
#SBATCH --mem 64G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step30a_eu15580_neighbours.%j.log
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 -u /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step30a_eu15580_neighbours.py
echo "exit $?"; date
