#!/bin/bash
#SBATCH -J step34f_cse2
#SBATCH -p stajichlab -c 8 --mem 48G --time=12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34f_cse2_loci.%j.log
# Step 34f: Cse2/CasB-like family membership and type I-E cas loci (see scripts/step34f_cse2_loci.py header).
# Tools: prodigal 2.6.3 (-p meta), samtools 1.22.1 (faidx on catalog/db/LsPS_AA_all.fasta.gz), minimap2 2.30
#   (-x asm10 --secondary=no -c, contigs without a supported-rep origin label only), hmmer 3.4 hmmsearch --cut_ga
#   vs Pfam-A 38.2 (db-pfam, staged to node-local scratch), MinCED 0.4.2 (condaenvs/minced, default settings).
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
D=$W/results/step34f_cse2
for f in $W/results/step34f_cse2_loci.tsv $W/results/fig_prep_cse2_locus_best.tsv $D; do
  if [ -e $f ]; then echo "REFUSING: $f exists"; exit 1; fi
done
hostname; date
$PY $W/scripts/step34f_cse2_loci.py members $D || { echo "exit 1 (members)"; exit 1; }
module load prodigal/2.6.3
prodigal -v 2>&1 | head -2
prodigal -p meta -i $D/contigs.fa -a $D/prodigal.faa -f gff -o $D/prodigal.gff -q || { echo "exit 1 (prodigal)"; exit 1; }
$PY $W/scripts/step34f_cse2_loci.py verify $D || exit 1
module load samtools/1.22.1
samtools faidx $W/catalog/db/LsPS_AA_all.fasta.gz -r $D/catalog_ids.txt > $D/catalog.faa || { echo "exit 1 (faidx)"; exit 1; }
$PY $W/scripts/step34f_cse2_loci.py compare $D || { echo "exit 1 (compare)"; exit 1; }
if [ -s $D/tomap.txt ]; then
  module load minimap2/2.30
  module load workspace/scratch
  samtools faidx $D/contigs.fa
  samtools faidx $D/contigs.fa -r $D/tomap.txt > $D/tomap.fa
  minimap2 -t 8 -x asm10 --secondary=no -c $W/catalog/fly/dm6.mmi $D/tomap.fa > $D/fly.paf
  cp $W/catalog/host/aRanSyl1.mmi $SCRATCH/host.mmi && minimap2 -t 8 -x asm10 --secondary=no -c $SCRATCH/host.mmi $D/tomap.fa > $D/frog.paf
  echo "mapped $(grep -c '>' $D/tomap.fa) contigs: fly rows $(wc -l < $D/fly.paf), frog rows $(wc -l < $D/frog.paf)"
fi
module load hmmer/3.4
module load db-pfam
echo "PFAM_DB $PFAM_DB"
module load workspace/scratch
mkdir -p $SCRATCH/pfam; for F in $PFAM_DB/Pfam-A.hmm*; do cp $F $SCRATCH/pfam/; done
hmmsearch --cut_ga --cpu 8 --tblout $D/hmmsearch.tblout --domtblout $D/hmmsearch.domtblout $SCRATCH/pfam/Pfam-A.hmm $D/prodigal.faa > /dev/null || { echo "exit 1 (hmmsearch)"; exit 1; }
/bigdata/stajichlab/lshad003/condaenvs/minced/bin/minced --version | head -1
/bigdata/stajichlab/lshad003/condaenvs/minced/bin/minced -gffFull $D/contigs.fa $D/minced.txt $D/minced.gff || { echo "exit 1 (minced)"; exit 1; }
$PY $W/scripts/step34f_cse2_loci.py loci $D
echo "exit $?"; date
