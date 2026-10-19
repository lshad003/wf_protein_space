#!/bin/bash
#SBATCH -J step34m_tank
#SBATCH -p stajichlab -c 2 --mem 48G --time=12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34m_tank_restriction.%j.log
# Step 34m wrapper: runs scripts/step34m_tank_restriction.py (tank restriction within clutch; WF22 one sample per animal)
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34m_tank_restriction.py
echo "exit $?"; date
