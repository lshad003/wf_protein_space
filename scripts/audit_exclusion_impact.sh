#!/bin/bash
#SBATCH -p stajichlab -c 2 --mem 8G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/audit_exclusion_impact.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/audit_exclusion_impact.py
