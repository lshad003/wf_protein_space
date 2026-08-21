#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 128gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step2e.%j.log
module load mmseqs2/13-45111
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
TMP=$SCRATCH/step2e_$SLURM_JOB_ID
mkdir -p $TMP
CPU=${SLURM_CPUS_ON_NODE:-32}

run () {  # name query target minid covmode
  mmseqs easy-search $2 $3 $WPS/results/step2e_$1.m8 $TMP/$1 \
    --min-seq-id $4 -c 0.8 --cov-mode $5 --max-accept 1 -e 1e-5 -s 4 --threads $CPU > /dev/null 2>&1
  Q=$(grep -c '^>' $2)
  H=$(cut -f1 $WPS/results/step2e_$1.m8 | sort -u | wc -l)
  $PY -c "print('%-28s %7d/%7d = %6.2f%%' % ('$1', $H, $Q, 100.0*$H/$Q))"
}

echo "==== cov-mode 2 (coverage of QUERY) ===="
run control_cov2 $WPS/catalog/wf22_control_200k.faa $V1/db/LsFMGC_AA_95_rep.fasta 0.95 2
run new95_cov2   $WPS/catalog/subsample_200k.faa    $V1/db/LsFMGC_AA_95_rep.fasta 0.95 2
run new50_cov2   $WPS/catalog/subsample_200k.faa    $V1/db/LsFMGC_AA_50_rep.fasta 0.50 2
echo
echo "==== cov-mode 0 (both sequences, strictest) ===="
run control_cov0 $WPS/catalog/wf22_control_200k.faa $V1/db/LsFMGC_AA_95_rep.fasta 0.95 0
run new95_cov0   $WPS/catalog/subsample_200k.faa    $V1/db/LsFMGC_AA_95_rep.fasta 0.95 0
rm -rf $TMP
