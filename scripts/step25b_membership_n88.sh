#!/bin/bash
#SBATCH -p stajichlab -c 2 --mem 32G --time 6:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step25b.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step25b_membership_n88.py
