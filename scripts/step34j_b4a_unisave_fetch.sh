#!/bin/bash
#SBATCH -p stajichlab
#SBATCH --constraint=milan
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -t 3-00:00:00
#SBATCH -J step34j_b4a_unisave
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_b4a_unisave.%j.log
# Step 34j B4a: UniProt secondary accessions and UniSave last versions for the b4 gap (scripts/step34j_b4a_unisave_fetch.py).
# Resumable: resubmitting skips accessions already recorded. Optional $1: extra accession list (passed to the .py).
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34j_b4a_unisave_fetch.py $1
rc=$?; echo "exit $rc"; date; exit $rc
