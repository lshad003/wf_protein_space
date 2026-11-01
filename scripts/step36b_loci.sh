#!/bin/bash
#SBATCH -J step36b_loci
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem 48G
#SBATCH -t 2:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step36b_loci.%j.log
# Step 36b A: summary rows and locus plots of the 8 several-variant structure-confident EU families (scripts/step36b_loci.py).
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step36b_loci.py
rc=$?
echo "exit $rc"; date
exit $rc
