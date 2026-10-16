#!/bin/bash
#SBATCH -J step33g_map
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33g_map.%A_%a.log
# Step 33g: map the step33c contig lists (results/step33c_lists/<stem>.tsv: sets A, B, C, F; unchanged) to
# D. melanogaster dm6 (catalog/fly/dm6.mmi, 144 Mb, read in place) with minimap2 2.30 -x asm10 --secondary=no -c.
# Array index N = line N of results/step33c_lists/stems.txt (89: the 88, then UHM586.41010 for set B).
# Contigs extracted from the MEGAHIT assemblies into $SCRATCH. Output: results/step33g_paf/<stem>.fly.paf; skips if it exists.
W=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results
N=${SLURM_ARRAY_TASK_ID:-$1}
S=$(sed -n ${N}p $W/results/step33c_lists/stems.txt)
if [ -z "$S" ]; then echo "no stem for index $N"; exit 1; fi
O=$W/results/step33g_paf; mkdir -p $O
echo "index $N stem $S"; hostname; date
if [ -e $O/$S.fly.paf ]; then echo "exists, skipping"; exit 0; fi
module load minimap2/2.30
module load workspace/scratch
T=$SCRATCH/step33g_$S; mkdir -p $T
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - $FEC/$S/scaffolds/${S}_R.fa.gz $W/results/step33c_lists/$S.tsv $T/q.fa <<'PYEOF' || { echo "ABORT: extraction failed"; exit 1; }
import sys, gzip
asm, lst, oq = sys.argv[1:]
want = {l.split("\t")[0] for l in open(lst)}
n = 0; w = False
with gzip.open(asm, "rt") as fh, open(oq, "w") as q:
    for l in fh:
        if l[0] == ">": w = l[1:].split()[0] in want; n += w
        if w: q.write(l)
print("listed %d, extracted %d" % (len(want), n))
if n != len(want): sys.exit("ABORT: %d listed contigs not found" % (len(want) - n))
PYEOF
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no -c $W/catalog/fly/dm6.mmi $T/q.fa > $T/fly.paf
rc=$?
echo "minimap2 exit $rc"
if [ $rc -ne 0 ]; then rm -rf $T; exit 1; fi
cp $T/fly.paf $O/$S.fly.paf
rc=$?
echo "fly PAF rows $(wc -l < $O/$S.fly.paf), contigs $(cut -f1 $O/$S.fly.paf | sort -u | wc -l)"
rm -rf $T
echo "exit $rc"; date
exit $rc
