#!/bin/bash
#SBATCH -J step32a2_tiara_41010
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step32a2_tiara_41010.%j.log
# Step 32a2: Tiara on the excluded control UHM586.41010 (outside the 88-array), same settings as step32a.
# Used in step32b only to label supported 95% reps whose own contig is in UHM586.41010; the sample stays out of
# all member, family, abundance and cohort counts. Output: results/tiara/UHM586.41010.tiara.tsv; skips if it exists.
W=/bigdata/stajichlab/lshad003/wf_protein_space
S=UHM586.41010
IN=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results/$S/scaffolds/${S}_R.fa.gz
OUT=$W/results/tiara/$S.tiara.tsv
hostname; date
if [ -e $OUT ]; then echo "exists, skipping: $OUT"; exit 0; fi
if [ ! -s $IN ]; then echo "ABORT: no contigs for $S"; exit 1; fi
/bigdata/stajichlab/lshad003/condaenvs/tiara/bin/tiara -i $IN -o $OUT -m 1000 -t $SLURM_CPUS_PER_TASK --probabilities
rc=$?
echo "exit $rc"; date
if [ $rc -eq 0 ]; then echo "class counts:"; tail -n +2 $OUT | cut -f2 | sort | uniq -c; fi
exit $rc
