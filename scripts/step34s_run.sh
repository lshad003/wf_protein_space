#!/bin/bash
#SBATCH -p epyc --constraint=milan
# Step 34s: read-count analyses on one count matrix (usage: sbatch -J step34s_<label>_<part> -c N --mem M -t T
#   -o logs/step34s_<label>_<part>.%j.log scripts/step34s_run.sh <label> <matrix> <part>).
# Outputs in results/step34s/<label>/ (the scripts refuse to overwrite). Labels: old (reproduction check on
#   results/count_matrix_primary_n88.tsv.gz), c_id95, c_id90, c_id95_q1 (results/step34r/matrix_<label>_n88.tsv.gz).
# part py:  step34b2_final_stats.py (read share, prevalence, three-cohort presence), step16v_prevalence_n88.py
#           (Fig. 2d), step34k_final.py Part D (read retention, Fig. 1d)
# part r:   step34e_counts.py, step34e_hostfactor.R (strict_final_EU, strict_final; PERMANOVA and limma-voom),
#           step34e_summary.py, step34l_m23_permanova.R (WF22 months 2-3), fig_prep_pcoa.R (Fig. 5a,b)
# part rsens: part r without fig_prep_pcoa.R (sensitivity matrices)
# Seeds, filters and permutation schemes are those of the scripts; nothing is changed here. Exits with the status of
#   the first failing command, else 0 (CLAUDE.md rule).
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
LAB=$1; MAT=$2; PART=$3
O=$W/results/step34s/$LAB
[ -s "$MAT" ] || { echo "ABORT: matrix $MAT missing"; exit 1; }
mkdir -p $O
echo "label $LAB matrix $MAT part $PART host $(hostname) $(date)"
run() { echo "== $* ($(date))"; "$@"; rc=$?; echo "== exit $rc"; if [ $rc -ne 0 ]; then echo "ABORT at: $*"; exit $rc; fi; }
if [ "$PART" = "py" ]; then
  run $PY $W/scripts/step34b2_final_stats.py $MAT $O
  run $PY $W/scripts/step16v_prevalence_n88.py $MAT $O
  run $PY $W/scripts/step34k_final.py $MAT $O
elif [ "$PART" = "r" ] || [ "$PART" = "rsens" ]; then
  module load R
  Rscript -e 'cat(R.version.string, "limma", as.character(packageVersion("limma")), "edgeR", as.character(packageVersion("edgeR")), "vegan", as.character(packageVersion("vegan")), "\n")'
  run $PY $W/scripts/step34e_counts.py $MAT $O
  run Rscript $W/scripts/step34e_hostfactor.R strict_final_EU $O/step34e_counts_strict_final_EU.tsv.gz $O
  run Rscript $W/scripts/step34e_hostfactor.R strict_final $O/step34e_counts_strict_final.tsv.gz $O
  run $PY $W/scripts/step34e_summary.py $O
  run Rscript $W/scripts/step34l_m23_permanova.R $O $O/step34l_m23_permanova.tsv
  if [ "$PART" = "r" ]; then run Rscript $W/scripts/fig_prep_pcoa.R $O $O/fig_pcoa_wf22.tsv; fi
else
  echo "ABORT: unknown part $PART"; exit 1
fi
echo "DONE $LAB $PART $(date)"
exit 0
