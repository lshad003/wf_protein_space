#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 2:00:00
#SBATCH -J step34t2
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34t2_tank_run_check.%j.log
# Step 34t2: tank vs sequencing run, within-run and exact p, saved nulls (scripts/step34t2_tank_run_check.py).
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34t2_tank_run_check.py
rc=$?; echo "exit $rc"; date; exit $rc
