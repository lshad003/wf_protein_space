#!/bin/bash
#SBATCH -p epyc -c 4 --mem 48G -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34s_stp1717_check.%j.log
# Step 34s check: runs scripts/step34s_stp1717_check.R (STP1717.1 LOAO and months 2-3 refits on c_id95).
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step34s_stp1717_check.R "$@"
rc=$?
echo "exit $rc $(date)"
exit $rc
