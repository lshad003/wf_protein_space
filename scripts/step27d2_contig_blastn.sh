#!/bin/bash
#SBATCH -p short -c 8 --mem 16G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step27d2.%j.log
# Step 27d2: blastn of the one bacterial-context contig (UHM611.41020 k141_15734, step27c) against the frozen
# Basidiobolus genomes. Contig from the MEGAHIT scaffolds the catalog genes were predicted from
# (Fecal/results/UHM611.41020/scaffolds/UHM611.41020_R.fa.gz, as in Fecal/pipeline/10_predict_by_assembl.sh).
# Prodigal (Fecal/Proteins/by_assembly/UHM611.41020.gff.gz): gene _1 1-315 (+, partial 10, neighbour, Bacteria),
# gene _2 339-533 (+, partial 01, the target). Contig length 533.
# Genomes: all *_softmasked.fasta.gz in the frozen set. The prompt names five; the set holds six (UHM516_7697 has a
# genome but no proteins file, so step27a did not search it). All six are searched and reported per genome.
# blastn defaults (megablast task not forced: -task blastn), -evalue 1e-5, softmasking ignored for seeding
# (-soft_masking false), all HSPs kept.
W=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
B=/bigdata/stajichlab/shared/projects/Herptile/Basidiobolus/FINAL_FROZEN_GENOMES_BASIDIOBOLUS_AUG18_2025
OUT=$W/results/basid_contig_k141_15734_blastn.tsv
if [ -e $OUT ]; then echo "ABORT: $OUT exists"; exit 1; fi
module load ncbi-blast/2.14.1+
module load workspace/scratch
T=$SCRATCH/step27d2; mkdir -p $T
zcat $FEC/results/UHM611.41020/scaffolds/UHM611.41020_R.fa.gz \
  | awk '/^>/{p=($1==">k141_15734")} p' | sed 's/^>k141_15734/>UHM611.41020__k141_15734/' > $T/contig.fa
echo "contig records: $(grep -c '^>' $T/contig.fa); length: $(grep -v '^>' $T/contig.fa | tr -d '\n' | wc -c)"
echo -e "genome\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > $OUT
for g in $B/*_FINAL/*_softmasked.fasta.gz; do
  n=$(basename $g _softmasked.fasta.gz)
  zcat $g > $T/$n.fa
  makeblastdb -in $T/$n.fa -dbtype nucl -out $T/$n > /dev/null
  blastn -task blastn -query $T/contig.fa -db $T/$n -evalue 1e-5 -soft_masking false -num_threads 8 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    | sed "s/^/$n\t/" >> $OUT
  echo "$n: $(grep -c "^$n	" $OUT) HSPs"
done
rm -rf $T
echo "DONE step27d2"
