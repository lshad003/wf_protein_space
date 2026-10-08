#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 64G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step21b_n88.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step21b_neighbours_n88.py
