#!/bin/bash
#SBATCH -J step29d_nb
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem 16G
#SBATCH -t 2:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step29d_fold4_neighbours.%j.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step29d_fold4_neighbours.py
