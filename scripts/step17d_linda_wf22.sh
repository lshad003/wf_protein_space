#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 120G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step17d.log
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step17d_linda_wf22.R
