#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 16G --time 2:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step19a_n88.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step19a_core_unknown_n88.py
