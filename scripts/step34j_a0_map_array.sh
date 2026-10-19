#!/bin/bash
#SBATCH -J step34j_map
#SBATCH -p stajichlab -c 8 --mem 48G -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_a0_map.%A_%a.log
# Step 34j A0: map control member contigs without a supported-rep label and without a Step 34g PAF
# (results/step34j/ctl_maplists/<stem>.txt) to dm6 and the wood frog, settings as scripts/step34g_map_array.sh:
# minimap2 2.30 -x asm10 --secondary=no -c; frog index staged to node-local scratch.
# Array index N = line N of results/step34j/ctl_stems_to_map.txt. Outputs results/step34j/ctl_paf/<stem>.{fly,frog}.paf;
# refuses to overwrite.
W=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results
D=$W/results/step34j
S=$(sed -n ${SLURM_ARRAY_TASK_ID}p $D/ctl_stems_to_map.txt)
[ -z "$S" ] && { echo "no stem"; exit 1; }
O=$D/ctl_paf; mkdir -p $O
echo "index $SLURM_ARRAY_TASK_ID stem $S"; hostname; date
if [ -e $O/$S.fly.paf ] || [ -e $O/$S.frog.paf ]; then echo "REFUSING: PAF for $S exists"; exit 1; fi
module load minimap2/2.30
module load workspace/scratch
T=$SCRATCH/step34j_map_$S; mkdir -p $T
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - $FEC/$S/scaffolds/${S}_R.fa.gz $D/ctl_maplists/$S.txt $T/q.fa <<'PYEOF' || { echo "ABORT: extraction"; exit 1; }
import sys, gzip
a, lst, oq = sys.argv[1:]
want = {l.strip() for l in open(lst) if l.strip()}; n = 0; w = False
with gzip.open(a, "rt") as fh, open(oq, "w") as q:
    for l in fh:
        if l[0] == ">": w = l[1:].split()[0] in want; n += w
        if w: q.write(l)
print("listed %d extracted %d" % (len(want), n))
if n != len(want): sys.exit(1)
PYEOF
cp $W/catalog/host/aRanSyl1.mmi $T/host.mmi || exit 1
minimap2 -t 8 -x asm10 --secondary=no -c $T/host.mmi $T/q.fa > $T/frog.paf; r1=$?
minimap2 -t 8 -x asm10 --secondary=no -c $W/catalog/fly/dm6.mmi $T/q.fa > $T/fly.paf; r2=$?
[ $r1 -ne 0 ] || [ $r2 -ne 0 ] && { echo "minimap2 failed $r1 $r2"; rm -rf $T; exit 1; }
cp $T/frog.paf $O/$S.frog.paf && cp $T/fly.paf $O/$S.fly.paf
echo "frog rows $(wc -l < $O/$S.frog.paf) fly rows $(wc -l < $O/$S.fly.paf)"
rm -rf $T; echo "exit 0"; date
