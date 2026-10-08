#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 48G --time 8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step20c_n88.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step20c_gene_tax_n88.py
