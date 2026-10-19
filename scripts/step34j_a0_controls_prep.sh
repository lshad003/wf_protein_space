#!/bin/bash
#SBATCH -J step34j_a0
#SBATCH -p epyc
#SBATCH -c 2
#SBATCH --mem 64G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_a0_controls_prep.%j.log
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34j_a0_controls_prep.py
echo "exit $?"; date
