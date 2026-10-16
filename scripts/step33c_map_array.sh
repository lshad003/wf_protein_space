#!/bin/bash
#SBATCH -J step33c_map
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 64G
#SBATCH -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33c_map.%A_%a.log
# Step 33c: map listed contigs (results/step33c_lists/<stem>.tsv, step33c_prep.py) to the wood frog genome
# (catalog/host/aRanSyl1.mmi, GCA_028564925.1) and the set-A contigs also to the six Basidiobolus genomes
# (catalog/basid/basid6.mmi). minimap2 2.30 -x asm10 --secondary=no. Array index N = line N of
# results/step33c_lists/stems.txt (89: the 88, then UHM586.41010). The host index (8.8 GB) is staged to node-local
# $SCRATCH. Contigs are extracted from the MEGAHIT assemblies (as step32a) into $SCRATCH.
# Outputs: results/step33c_paf/<stem>.host.paf and <stem>.basid.paf; skips the task if both exist.
W=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
N=${SLURM_ARRAY_TASK_ID:-$1}
S=$(sed -n ${N}p $W/results/step33c_lists/stems.txt)
if [ -z "$S" ]; then echo "no stem for index $N"; exit 1; fi
O=$W/results/step33c_paf; mkdir -p $O
echo "index $N stem $S"; hostname; date
if [ -e $O/$S.host.paf ] && [ -e $O/$S.basid.paf ]; then echo "exists, skipping"; exit 0; fi
module load minimap2/2.30
module load workspace/scratch
T=$SCRATCH/step33c_$S; mkdir -p $T
IN=$FEC/$S/scaffolds/${S}_R.fa.gz
if [ ! -s $IN ]; then echo "ABORT: no contigs for $S"; exit 1; fi
$PY - $IN $W/results/step33c_lists/$S.tsv $T/q.fa $T/qA.fa <<'PYEOF' || { echo "ABORT: extraction failed"; exit 1; }
import sys, gzip
asm, lst, oq, oa = sys.argv[1:]
want = dict(l.rstrip("\n").split("\t") for l in open(lst))
n = na = 0; w = None
with gzip.open(asm, "rt") as fh, open(oq, "w") as q, open(oa, "w") as a:
    for l in fh:
        if l[0] == ">":
            c = l[1:].split()[0]; w = want.get(c)
            if w is not None: n += 1; na += "A" in w
        if w is not None:
            q.write(l)
            if "A" in w: a.write(l)
print("listed %d, extracted %d, of them set A %d" % (len(want), n, na))
if n != len(want): sys.exit("ABORT: %d listed contigs not found" % (len(want) - n))
PYEOF
cp $W/catalog/host/aRanSyl1.mmi $T/host.mmi || { echo "ABORT: staging failed"; exit 1; }
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no $T/host.mmi $T/q.fa > $T/host.paf
rc1=$?
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no $W/catalog/basid/basid6.mmi $T/qA.fa > $T/basid.paf
rc2=$?
echo "minimap2 exit host $rc1 basid $rc2"
if [ $rc1 -ne 0 ] || [ $rc2 -ne 0 ]; then rm -rf $T; exit 1; fi
cp $T/host.paf $O/$S.host.paf && cp $T/basid.paf $O/$S.basid.paf
rc=$?
echo "host PAF rows $(wc -l < $O/$S.host.paf), contigs $(cut -f1 $O/$S.host.paf | sort -u | wc -l); basid PAF rows $(wc -l < $O/$S.basid.paf)"
rm -rf $T
echo "exit $rc"; date
exit $rc
