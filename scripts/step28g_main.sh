#!/bin/bash
#SBATCH -p epyc -c 32 --mem 64G -t 24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28g.%j.log
# Step 28g: ProstT5 (CPU) + Foldseek search of the 15,580 recurrent EU-only family reps, as step28e
# (search vs PDB and CATH50, -e 1e-1); convertalis fields as step28e2 (no prob, TM or lddt).
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
M=/bigdata/stajichlab/lshad003/software/foldseek_gpu/weights/prostt5
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
W=/bigdata/stajichlab/lshad003/wf_protein_space
O=$W/results/prostt5_main; mkdir -p $O/tmp
if [ -e $O/eu15580_vs_PDB.tsv ]; then echo "ABORT: outputs exist"; exit 1; fi
date
$F createdb $W/results/prostt5_main_eu15580.faa $O/eu15580 --prostt5-model $M --threads 32
date
for t in PDB CATH50; do
  $F search $O/eu15580 $DB/$t $O/res_$t $O/tmp -e 1e-1 --threads 32
  $F convertalis $O/eu15580 $DB/$t $O/res_$t $O/eu15580_vs_$t.tsv \
    --format-output "query,target,evalue,bits,fident,alnlen,qlen,tlen,qcov,tcov"
  echo "$t: hit rows $(wc -l < $O/eu15580_vs_$t.tsv); queries with a hit $(cut -f1 $O/eu15580_vs_$t.tsv | sort -u | wc -l)"
done
date
echo "DONE step28g"
