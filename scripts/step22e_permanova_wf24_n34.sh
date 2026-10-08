#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 48G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step22e.log
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step22e_permanova_wf24_n34.R
