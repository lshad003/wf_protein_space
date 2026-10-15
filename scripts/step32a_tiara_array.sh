#!/bin/bash
#SBATCH -J step32a_tiara
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step32a_tiara.%A_%a.log
# Step 32a: Tiara prokaryote/eukaryote classification of the assembled contigs of the 88 metagenomes.
# Stems: data rows of metadata/samples_88_treatment.tsv (column 1); array index N = data row N (1-88).
# Contigs: as step1b/step1c, $FEC/<stem>/scaffolds/<stem>_R.fa(.gz) (all 88 are .fa.gz).
# tiara (env /bigdata/stajichlab/lshad003/condaenvs/tiara) -m 1000 -t 8 --probabilities.
# Output: results/tiara/<stem>.tiara.tsv; skips the stem if it exists. Run one stem without Slurm: bash <script> N
W=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results
TIARA=/bigdata/stajichlab/lshad003/condaenvs/tiara/bin/tiara
N=${SLURM_ARRAY_TASK_ID:-$1}
if [ -z "$N" ]; then echo "need array index"; exit 1; fi
S=$(tail -n +2 $W/metadata/samples_88_treatment.tsv | sed -n ${N}p | cut -f1)
if [ -z "$S" ]; then echo "no stem for index $N"; exit 1; fi
OUT=$W/results/tiara/$S.tiara.tsv
echo "index $N stem $S"; hostname; date
if [ -e $OUT ]; then echo "exists, skipping: $OUT"; exit 0; fi
IN=$FEC/$S/scaffolds/${S}_R.fa
if [ ! -s $IN ]; then IN=$IN.gz; fi
if [ ! -s $IN ]; then echo "ABORT: no contigs for $S"; exit 1; fi
mkdir -p $W/results/tiara
$TIARA -i $IN -o $OUT -m 1000 -t ${SLURM_CPUS_PER_TASK:-8} --probabilities
rc=$?
echo "exit $rc"; date
if [ $rc -eq 0 ]; then echo "class counts:"; tail -n +2 $OUT | cut -f2 | sort | uniq -c; fi
exit $rc
