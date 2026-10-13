#!/bin/bash
#SBATCH -p epyc -c 32 --mem 64G --time 24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28e.%j.log
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
M=/bigdata/stajichlab/lshad003/software/foldseek_gpu/weights/prostt5
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
W=/bigdata/stajichlab/lshad003/wf_protein_space
O=$W/results/prostt5_pilot; mkdir -p $O/tmp
if [ -e $O/eu740_vs_PDB.tsv ]; then echo "ABORT: outputs exist"; exit 1; fi
date
$F createdb $W/results/prostt5_pilot_eu740.faa $O/eu740 --prostt5-model $M --threads 32
date
for t in PDB CATH50; do
  $F search $O/eu740 $DB/$t $O/res_$t $O/tmp -e 1e-1 --threads 32
  $F convertalis $O/eu740 $DB/$t $O/res_$t $O/eu740_vs_$t.tsv \
    --format-output "query,target,evalue,bits,prob,fident,alnlen,qlen,tlen,qcov,tcov"
  echo "$t: hit rows $(wc -l < $O/eu740_vs_$t.tsv); queries with a hit $(cut -f1 $O/eu740_vs_$t.tsv | sort -u | wc -l)"
done
date
echo "DONE step28e"
