#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 16G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16u_n88.log
R=/bigdata/stajichlab/lshad003/wf_protein_space/results
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step16u_matrix_n88.py || exit 1
cmp <(zcat $R/count_matrix_primary_n88.check.tsv.gz) <(zcat $R/count_matrix_primary_n88.tsv.gz) \
  && echo "primary n88 re-derived matrix identical to step25a output" && rm $R/count_matrix_primary_n88.check.tsv.gz
