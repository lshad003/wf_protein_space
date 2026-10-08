#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 64G --time 8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step24b_n88.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step24b_accum_n88.py
