#!/bin/bash
#SBATCH -p epyc --constraint=milan -c 8 --mem 160G -t 06:00:00
#SBATCH -J step34s_month_sharedruns -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34s_month_sharedruns.%j.log
# Step 34s sensitivity: WF22 month PERMANOVA within the two shared sequencing runs (resources and module as step34s_run.sh part r).
W=/bigdata/stajichlab/lshad003/wf_protein_space
echo "host $(hostname) $(date)"
module load R
Rscript -e 'cat(R.version.string, "edgeR", as.character(packageVersion("edgeR")), "vegan", as.character(packageVersion("vegan")), "permute", as.character(packageVersion("permute")), "\n")'
Rscript $W/scripts/step34s_month_sharedruns.R
rc=$?
echo "exit $rc $(date)"
exit $rc
