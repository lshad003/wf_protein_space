#!/bin/bash
#SBATCH -J step33b_host_index
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 96G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33b_host_index.%j.log
# Step 33b: minimap2 2.30 index of the wood frog genome (Lithobates sylvaticus GCA_028564925.1, aRanSyl1.merge,
# catalog/host/, see catalog/host/SOURCE.txt) for mapping assembled contigs to the host.
# Preset -x asm10 (k=19, w=19; the same index serves asm5 and asm10 mapping). Short-read mapping needs a separate
# -x sr index. -I 8G keeps the 5.15 Gb genome in a single index part (a multi-part index degrades mapping).
# Output: catalog/host/aRanSyl1.mmi; refuses to overwrite.
W=/bigdata/stajichlab/lshad003/wf_protein_space
H=$W/catalog/host
IN=$H/GCA_028564925.1_aRanSyl1.merge_genomic.fna.gz
OUT=$H/aRanSyl1.mmi
if [ -e $OUT ]; then echo "REFUSING: $OUT exists"; exit 1; fi
if [ ! -s $IN ]; then echo "ABORT: missing $IN"; exit 1; fi
module load minimap2/2.30
hostname; date; minimap2 --version
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 -I 8G -d $OUT $IN
rc=$?
echo "exit $rc"; date; ls -la $OUT
exit $rc
