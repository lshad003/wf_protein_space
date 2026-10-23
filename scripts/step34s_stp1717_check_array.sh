#!/bin/bash
#SBATCH -p epyc -c 4 --mem 48G -t 8:00:00 --array=0-11
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34s_stp1717_check_array.%A_%a.log
# Step 34s check: one refit per task of scripts/step34s_stp1717_check_array.R (full, 10 leave-one-animal-out, months 2-3).
module load R
mkdir -p /bigdata/stajichlab/lshad003/wf_protein_space/results/step34s/stp1717_check/loao
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34s_stp1717_check_array.R
rc=$?
echo "exit $rc $(date)"
exit $rc
