#!/bin/bash
#SBATCH -p short -w r02 -c 4 --mem 16G -t 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28e2.%j.log
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
O=/bigdata/stajichlab/lshad003/wf_protein_space/results/prostt5_pilot
for t in PDB CATH50; do
  ls $O/res_$t.index >/dev/null || { echo "MISSING res_$t"; continue; }
  $F convertalis $O/eu740 $DB/$t $O/res_$t $O/eu740_vs_$t.tsv \
    --format-output "query,target,evalue,bits,fident,alnlen,qlen,tlen,qcov,tcov"
  echo "$t: hit rows $(wc -l < $O/eu740_vs_$t.tsv); queries with a hit $(cut -f1 $O/eu740_vs_$t.tsv | sort -u | wc -l)"
done
echo "DONE step28e2"
