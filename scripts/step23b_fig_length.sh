#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 24G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step23b.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step23b_fig_length.py
