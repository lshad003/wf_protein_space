#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 2:00:00
#SBATCH -J step34t3
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34t3_tank_audit.%j.log
# Step 34t3: WF24 EM4 tank audit (scripts/step34t3_tank_audit.py).
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34t3_tank_audit.py
rc=$?; echo "exit $rc"; date; exit $rc
