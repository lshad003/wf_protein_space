#!/bin/bash
#SBATCH -J step28h_nostar
#SBATCH -p epyc
#SBATCH -c 32
#SBATCH --mem 64G
#SBATCH -t 2-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28h_nostar.%j.log
# Step 28h: rerun of the Step 28 ProstT5 (CPU) + Foldseek searches without the terminal '*' (step28h_inputs.py).
# Same binary, ProstT5 weights, databases and settings as step28e/28f/28g: createdb --prostt5-model, search vs PDB
# and CATH50 -e 1e-1, convertalis with the step28e2/28f/28g fields. Outputs in results/prostt5_{main,pilot,calib}_nostar;
# refuses to run if any of them exists. Then scripts/step28h_compare.py; exits with its exit code.
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
M=/bigdata/stajichlab/lshad003/software/foldseek_gpu/weights/prostt5
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
R=$W/results
for d in $R/prostt5_main_nostar $R/prostt5_pilot_nostar $R/prostt5_calib_nostar; do
  if [ -e $d ]; then echo "REFUSING: $d exists"; exit 1; fi
done
# query FASTA : output dir : db prefix
SETS="$R/prostt5_main_eu15580_nostar.faa:$R/prostt5_main_nostar:eu15580
$R/prostt5_pilot_eu740_nostar.faa:$R/prostt5_pilot_nostar:eu740
$R/prostt5_calib_known_lenmatched_nostar.faa:$R/prostt5_calib_nostar:known740
$R/prostt5_decoy_eu740rev_nostar.faa:$R/prostt5_calib_nostar:decoy740"
hostname; date; $F version
for s in $SETS; do
  IFS=: read faa O q <<< "$s"
  mkdir -p $O/tmp
  echo "== $q ($(grep -c '^>' $faa) sequences)"; date
  $F createdb $faa $O/$q --prostt5-model $M --threads $SLURM_CPUS_PER_TASK
  for t in PDB CATH50; do
    $F search $O/$q $DB/$t $O/${q}_res_$t $O/tmp -e 1e-1 --threads $SLURM_CPUS_PER_TASK
    $F convertalis $O/$q $DB/$t $O/${q}_res_$t $O/${q}_vs_$t.tsv \
      --format-output "query,target,evalue,bits,fident,alnlen,qlen,tlen,qcov,tcov"
    if [ ! -s $O/${q}_vs_$t.tsv ]; then echo "ABORT: no hits written for $q vs $t"; exit 1; fi
    echo "$q $t: hit rows $(wc -l < $O/${q}_vs_$t.tsv); queries with a hit $(cut -f1 $O/${q}_vs_$t.tsv | sort -u | wc -l)"
  done
done
date
$PY $W/scripts/step28h_compare.py
rc=$?
echo "exit $rc"; date
exit $rc
