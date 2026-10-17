#!/bin/bash
#SBATCH -J step34e_host
#SBATCH -p stajichlab -c 8 --mem 160G --time=2-00:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34e_hostfactor.%j.log
# Step 34e: host-factor tests on strict_final and strict_final EU (designs of Steps 16 and 20; see step34e_hostfactor.R).
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
hostname; date
module load R
Rscript -e 'cat("limma", as.character(packageVersion("limma")), "edgeR", as.character(packageVersion("edgeR")), "vegan", as.character(packageVersion("vegan")), "\n")'
$PY $W/scripts/step34e_counts.py || { echo "exit 1 (counts)"; exit 1; }
Rscript $W/scripts/step34e_hostfactor.R strict_final_EU $W/results/step34e_counts_strict_final_EU.tsv.gz || { echo "exit 1 (EU)"; exit 1; }
Rscript $W/scripts/step34e_hostfactor.R strict_final $W/results/step34e_counts_strict_final.tsv.gz || { echo "exit 1 (strict_final)"; exit 1; }
$PY $W/scripts/step34e_summary.py
echo "exit $?"; date
