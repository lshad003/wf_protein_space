#!/bin/bash
#SBATCH -J step34a_map
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 64G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34a_map.%A_%a.log
# Step 34a: map the supported-rep contigs not mapped before (results/step34a_lists/<stem>.txt, step34a_prep.py) to the
# wood frog (catalog/host/aRanSyl1.mmi, staged to node-local $SCRATCH) and D. melanogaster (catalog/fly/dm6.mmi)
# indexes with the step33c2/step33g settings: minimap2 2.30 -x asm10 --secondary=no -c. Array index N = line N of
# results/step34a_lists/stems.txt (89). Contigs extracted from the MEGAHIT assemblies into $SCRATCH.
# Outputs: results/step34a_paf/<stem>.frog.paf and <stem>.fly.paf; skips the task if both exist.
W=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results
N=${SLURM_ARRAY_TASK_ID:-$1}
S=$(sed -n ${N}p $W/results/step34a_lists/stems.txt)
if [ -z "$S" ]; then echo "no stem for index $N"; exit 1; fi
O=$W/results/step34a_paf; mkdir -p $O
echo "index $N stem $S"; hostname; date
if [ -e $O/$S.frog.paf ] && [ -e $O/$S.fly.paf ]; then echo "exists, skipping"; exit 0; fi
module load minimap2/2.30
module load workspace/scratch
T=$SCRATCH/step34a_$S; mkdir -p $T
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - $FEC/$S/scaffolds/${S}_R.fa.gz $W/results/step34a_lists/$S.txt $T/q.fa <<'PYEOF' || { echo "ABORT: extraction failed"; exit 1; }
import sys, gzip
asm, lst, oq = sys.argv[1:]
want = {l.strip() for l in open(lst) if l.strip()}
n = 0; w = False
with gzip.open(asm, "rt") as fh, open(oq, "w") as q:
    for l in fh:
        if l[0] == ">": w = l[1:].split()[0] in want; n += w
        if w: q.write(l)
print("listed %d, extracted %d" % (len(want), n))
if n != len(want): sys.exit("ABORT: %d listed contigs not found" % (len(want) - n))
PYEOF
if [ ! -s $T/q.fa ]; then echo "no contigs to map"; : > $O/$S.frog.paf; : > $O/$S.fly.paf; rm -rf $T; echo "exit 0"; exit 0; fi
cp $W/catalog/host/aRanSyl1.mmi $T/host.mmi || { echo "ABORT: staging failed"; exit 1; }
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no -c $T/host.mmi $T/q.fa > $T/frog.paf
rc1=$?
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no -c $W/catalog/fly/dm6.mmi $T/q.fa > $T/fly.paf
rc2=$?
echo "minimap2 exit frog $rc1 fly $rc2"
if [ $rc1 -ne 0 ] || [ $rc2 -ne 0 ]; then rm -rf $T; exit 1; fi
cp $T/frog.paf $O/$S.frog.paf && cp $T/fly.paf $O/$S.fly.paf
rc=$?
echo "frog PAF rows $(wc -l < $O/$S.frog.paf); fly PAF rows $(wc -l < $O/$S.fly.paf)"
rm -rf $T
echo "exit $rc"; date
exit $rc
