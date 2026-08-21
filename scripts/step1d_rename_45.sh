#!/bin/bash
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
OLD=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly
NEW=$WPS/catalog/by_assembly_new
PEP=$WPS/catalog/input_pep
CDS=$WPS/catalog/input_cds
mkdir -p $PEP $CDS

cat $WPS/metadata/wf23_stems.txt $WPS/metadata/wf24_stems.txt > $WPS/metadata/new45_stems.txt
echo "stems to process (expect 45): $(wc -l < $WPS/metadata/new45_stems.txt)"

while read S; do
  if [ -s "$NEW/$S.aa.fa.gz" ]; then SRC=$NEW
  elif [ -s "$OLD/$S.aa.fa.gz" ]; then SRC=$OLD
  else echo "NOSOURCE $S"; continue; fi
  export BIOSAMPLE=$S
  pigz -dc $SRC/$S.aa.fa.gz  | perl -p -e 's/^>(\S+).+/>$ENV{BIOSAMPLE}__$1/' > $PEP/$S.fasta
  pigz -dc $SRC/$S.cds.fa.gz | perl -p -e 's/^>(\S+).+/>$ENV{BIOSAMPLE}__$1/' > $CDS/$S.fasta
  echo "done $S from $(basename $SRC)"
done < $WPS/metadata/new45_stems.txt
