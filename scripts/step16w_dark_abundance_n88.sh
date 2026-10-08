#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 16G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16w_n88.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step16w_dark_abundance_n88.py
