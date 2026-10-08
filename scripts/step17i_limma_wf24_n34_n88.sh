#!/bin/bash
#SBATCH -p stajichlab -c 8 --mem 64G --time 8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step17i_n88.log
module load R
Rscript /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step17i_limma_wf24_n34_n88.R
