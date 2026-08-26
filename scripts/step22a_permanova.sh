#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 48G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step22a.log
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step22a_permanova.R
