#!/bin/bash
#SBATCH -J step36_context
#SBATCH -p epyc
#SBATCH -c 2
#SBATCH --mem 96G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step36_context.%j.log
# Step 36: genomic context of the 4,509 recurrent EU families and the 4,509 K controls (scripts/step36_context.py).
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step36_context.py
rc=$?
echo "exit $rc"; date
exit $rc
