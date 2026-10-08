#!/bin/bash
#SBATCH -p stajichlab -c 8 --mem 64G --time 8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step17g_n88.log
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step17g_limma_wf22_n88.R
