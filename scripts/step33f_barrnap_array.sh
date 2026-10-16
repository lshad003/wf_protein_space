#!/bin/bash
#SBATCH -J step33f_barrnap
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33f_barrnap.%A_%a.log
# Step 33f: rRNA prediction on the Tiara-eukaryotic contigs >= 1 kb (eukarya + organelle, results/tiara/<stem>.tiara.tsv)
# of the 88 samples. barrnap 0.9 --kingdom euk (defaults: --evalue 1e-06, --reject 0.25, --lencutoff 0.8).
# Array index N = data row N of metadata/samples_88_treatment.tsv. Contigs extracted from the MEGAHIT assembly
# (as step32a) into $SCRATCH. Outputs: results/step33f_barrnap/<stem>.gff and <stem>.rrna.fa (all rRNA hits:
# 18S, 28S, 5.8S, 5S); skips the stem if both exist.
W=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results
N=${SLURM_ARRAY_TASK_ID:-$1}
S=$(tail -n +2 $W/metadata/samples_88_treatment.tsv | sed -n ${N}p | cut -f1)
if [ -z "$S" ]; then echo "no stem for index $N"; exit 1; fi
O=$W/results/step33f_barrnap; mkdir -p $O
echo "index $N stem $S"; hostname; date
if [ -e $O/$S.gff ] && [ -e $O/$S.rrna.fa ]; then echo "exists, skipping"; exit 0; fi
# barrnap module (barrnap/0.9) only loads miniconda3 and activates an env, which fails in non-login batch shells
# (first array 29988288: 'barrnap: command not found'); the env's bin (barrnap, its perl, nhmmer, bedtools) is put on PATH.
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin:$PATH
barrnap --version 2>&1 || { echo "ABORT: barrnap not runnable"; exit 1; }
module load workspace/scratch
T=$SCRATCH/step33f_$S; mkdir -p $T
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - $FEC/$S/scaffolds/${S}_R.fa.gz $W/results/tiara/$S.tiara.tsv $T/euk.fa <<'PYEOF' || { echo "ABORT: extraction failed"; exit 1; }
import sys, gzip
asm, tia, out = sys.argv[1:]
want = set()
with open(tia) as fh:
    next(fh)
    for l in fh:
        f = l.split("\t", 2)
        if f[1] in ("eukarya", "organelle"): want.add(f[0].split()[0])
n = 0; w = False
with gzip.open(asm, "rt") as fh, open(out, "w") as o:
    for l in fh:
        if l[0] == ">": w = l[1:].split()[0] in want; n += w
        if w: o.write(l)
print("Tiara-eukaryotic contigs >= 1 kb: %d, extracted %d" % (len(want), n))
if n != len(want): sys.exit("ABORT: %d contigs not found" % (len(want) - n))
PYEOF
if [ ! -s $T/euk.fa ]; then echo "no Tiara-eukaryotic contigs"; : > $O/$S.gff; : > $O/$S.rrna.fa; rm -rf $T; exit 0; fi
barrnap --kingdom euk --threads $SLURM_CPUS_PER_TASK --outseq $T/rrna.fa $T/euk.fa > $T/rrna.gff
rc=$?
echo "barrnap exit $rc"
if [ $rc -ne 0 ]; then rm -rf $T; exit 1; fi
cp $T/rrna.gff $O/$S.gff && cp $T/rrna.fa $O/$S.rrna.fa
rc=$?
echo "rRNA hits: $(grep -vc '^#' $O/$S.gff) | 18S: $(grep -c 'Name=18S_rRNA' $O/$S.gff)"
rm -rf $T
echo "exit $rc"; date
exit $rc
