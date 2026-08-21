#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 96 --mem 384gb --time=7-00:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step2_cluster.%j.log
module load mmseqs2/13-45111
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
V1PEP=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1/input_pep
MEM=350G
CPU=${SLURM_CPUS_ON_NODE:-96}

# EDIT THIS LINE: wf22_stems.txt (50) or wf22_stems_44.txt (44)
WF22SET=$WPS/metadata/wf22_stems_44.txt

mkdir -p $WPS/catalog/db
AADB=$WPS/catalog/db/LsPS_AA_all.fasta

if [ ! -s $AADB.gz ]; then
  echo "assembling input from $(wc -l < $WF22SET) WF22 + 45 new"
  for S in $(cat $WF22SET); do cat $V1PEP/$S.fasta; done  > $AADB
  cat $WPS/catalog/input_pep/*.fasta                     >> $AADB
  echo "input proteins: $(grep -c '^>' $AADB)"
  echo "input metagenomes: $(grep '^>' $AADB | sed 's/>//; s/__.*//' | sort -u | wc -l)"
  if grep -q -E "^>UHM739|^>UHM740" $AADB; then
     echo "ERROR: non-wood-frog sample in input"; exit 1
  fi
  module load samtools
  bgzip --threads $CPU $AADB
  samtools faidx $AADB.gz
fi

DB=$WPS/catalog/db/LsPS_AA
NEXP=$(wc -l < $AADB.gz.fai)
NGOT=0
[ -s $DB.lookup ] && NGOT=$(wc -l < $DB.lookup)
if [ "$NGOT" != "$NEXP" ]; then
  echo "createdb: have $NGOT of $NEXP, rebuilding"
  rm -f $DB $DB.dbtype $DB.index $DB.lookup $DB.source ${DB}_h ${DB}_h.dbtype ${DB}_h.index
  mmseqs createdb $AADB.gz $DB --compressed 1
  NGOT=$(wc -l < $DB.lookup)
fi
if [ "$NGOT" != "$NEXP" ]; then
  echo "ERROR: createdb incomplete, $NGOT of $NEXP"; exit 1
fi
echo "database verified: $NGOT sequences" 

for identity in 0.95 0.9 0.5; do
  NAME=$(perl -e "print $identity * 100")
  echo "=== clustering at $identity ($(date)) ==="
  mmseqs cluster $DB ${DB}_${NAME}_cluster $SCRATCH \
    --min-seq-id $identity -c 0.8 --cov-mode 1 \
    --threads $CPU --split-memory-limit $MEM --kmer-per-seq 80
  [ -f ${DB}_${NAME}_cluster_rep ] || \
    mmseqs createsubdb ${DB}_${NAME}_cluster $DB ${DB}_${NAME}_cluster_rep
  [ -s ${DB}_${NAME}_rep.fasta ] || \
    mmseqs convert2fasta ${DB}_${NAME}_cluster_rep ${DB}_${NAME}_rep.fasta
  echo "tier $NAME representatives: $(grep -c '^>' ${DB}_${NAME}_rep.fasta)"
done
echo "done $(date)"
