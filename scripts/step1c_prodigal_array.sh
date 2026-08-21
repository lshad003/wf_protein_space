#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 1 --mem 16gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/prodigal.%a.log
module load prodigal
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
INPUT=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results
OUTPUT=$WPS/catalog/by_assembly_new
SAMPFILE=$WPS/metadata/missing14_samples.csv
TMPD=$WPS/tmp
mkdir -p $OUTPUT $TMPD

N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then N=$1; fi
if [ -z $N ]; then echo "need array index"; exit 1; fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN
do
  BASE=$OUTPUT/${STRAIN}
  IN=$INPUT/$STRAIN/scaffolds/${STRAIN}_R.fa
  echo "running on $STRAIN base=$BASE infile=$IN"
  if [ -f $BASE.gff.gz ]; then echo "already done, skipping"; continue; fi
  if [[ ! -f $IN && -f $IN.gz ]]; then
     pigz -dc $IN.gz > $TMPD/${STRAIN}_R.fa
     IN=$TMPD/${STRAIN}_R.fa
  elif [[ ! -f $IN ]]; then
     echo "no assembly for $IN or $IN.gz"; continue
  fi
  prodigal -i $IN -f gff -o $BASE.gff -a $BASE.aa.fa -d $BASE.cds.fa -p meta
  pigz $BASE.gff $BASE.aa.fa $BASE.cds.fa
  rm -f $TMPD/${STRAIN}_R.fa
  echo "finished $STRAIN"
done
