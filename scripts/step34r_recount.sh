#!/bin/bash
#SBATCH -J step34r_recount
#SBATCH -p stajichlab -c 2 --mem 96G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34r_recount.%j.log
# Step 34r part 2: matrices (step34r_matrices.py) then the recount (step34r_recount.py). Exits with the failing status.
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
S=/bigdata/stajichlab/lshad003/wf_protein_space/scripts
date
$PY $S/step34r_matrices.py
rc=$?
if [ $rc -ne 0 ]; then echo "ABORT: step34r_matrices.py exit $rc"; exit $rc; fi
$PY $S/step34r_recount.py
rc=$?
echo "step34r_recount.py exit $rc"; date
exit $rc
