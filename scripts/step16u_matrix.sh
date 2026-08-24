#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 64G --time 8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16u.log
P=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
S=/bigdata/stajichlab/lshad003/wf_protein_space/scripts/step16u_matrix.py
$P $S 1 primary
$P $S 2 mapq10
echo DONE step16u
