#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH -t 2:00:00
#SBATCH -J step34t
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34t_eggmass_tank.%j.log
# Step 34t: WF22 egg-mass breadth and tank restriction at c_id95 (scripts/step34t_eggmass_tank.py).
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34t_eggmass_tank.py
rc=$?; echo "exit $rc"; date; exit $rc
