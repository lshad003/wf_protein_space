#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 32G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step12e.log
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
[ -s $WPS/results/known_union_ids.txt ] && [ -s $WPS/results/eggnog_hit_ids.txt ] \
  || { echo "GUARD FAIL: union id lists missing"; exit 1; }
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  $WPS/scripts/step12e_family_union.py
