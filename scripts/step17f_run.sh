#!/bin/bash
#SBATCH -p stajichlab -c 16 --mem 64G --time 24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step17f_%x.log
#SBATCH -J linda
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step17f_linda_one.R $TAG $MOD
