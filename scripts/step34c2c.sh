#!/bin/bash
#SBATCH -J step34c2c
#SBATCH -p stajichlab -c 2 --mem 96G --time=12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34c2c.%j.log
# Step 34c2c: part 1 (taxonomy domain parser fix, downstream summaries) then part 2 (strict_v2).
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
hostname; date
$PY /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34c2c_taxv2.py || { echo "exit 1 (part 1)"; exit 1; }
$PY /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34c2c_strict_v2.py
echo "exit $?"; date
