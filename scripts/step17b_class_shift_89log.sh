#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 4G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step17b_89.log
# Runs the unchanged scripts/step17b_class_shift.py to record its 89-sample output
# in a log; no earlier log exists. Writes no result file.
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step17b_class_shift.py
