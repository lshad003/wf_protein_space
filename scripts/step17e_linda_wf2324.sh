#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 48G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step17e.log
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step17e_linda_wf2324.R
