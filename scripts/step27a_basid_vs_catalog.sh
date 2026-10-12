#!/bin/bash
#SBATCH -p epyc -c 32 --mem 128G --time 24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step27a.%j.log
W=/bigdata/stajichlab/lshad003/wf_protein_space
B=/bigdata/stajichlab/shared/projects/Herptile/Basidiobolus/FINAL_FROZEN_GENOMES_BASIDIOBOLUS_AUG18_2025
OUT=$W/results/basid_vs_catalog_hits.tsv
if [ -e $OUT ]; then echo "ABORT: $OUT exists"; exit 1; fi
module load mmseqs2/13-45111
module load workspace/scratch
T=$SCRATCH/step27a; mkdir -p $T
rm -f $T/q.faa
for s in STP1710:Basidiobolus_STP1710 STP1717:Basidiobolus_STP1717_1 UHM207:Basidiobolus_UHM207_4505 \
         UHM260:Basidiobolus_UHM260_5136 UHM520:Basidiobolus_UHM520_7734; do
  d=${s%%:*}; f=${s#*:}
  sed "s/^>/>${d}|/" $B/${d}_UHMproject_genome_FINAL/$f.proteins.fa >> $T/q.faa
done
echo "queries: $(grep -c '^>' $T/q.faa)"
mmseqs createdb $T/q.faa $T/q
mmseqs createdb $W/catalog/db/LsPS_AA_95_rep.fasta $T/cat
mmseqs search $T/q $T/cat $T/res $T/tmp -e 1e-5 --min-seq-id 0.9 -c 0.8 --cov-mode 1 --threads 32
mmseqs convertalis $T/q $T/cat $T/res $OUT \
  --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov,evalue,bits"
echo "hit rows: $(wc -l < $OUT); queries with a hit: $(cut -f1 $OUT | sort -u | wc -l)"
echo "DONE step27a"
