#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 8G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step13h_89.log
# Runs the unchanged scripts/step13h_gu_eu_split.sh to record its 89-sample output
# in a log; the original was run on the login node and left no log.
bash /bigdata/stajichlab/lshad003/wf_protein_space/scripts/step13h_gu_eu_split.sh
