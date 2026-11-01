#!/bin/bash
#SBATCH -J step36b_hmm
#SBATCH -p epyc
#SBATCH -c 16
#SBATCH --mem 24G
#SBATCH -t 7-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step36b_hmm.%A_%a.log
# Step 36b: Pfam near-miss scan. HMMER 3.4 hmmsearch against Pfam-A 38.2 (db-pfam module, the catalog annotation
# release; staged to node-local scratch), --max (all acceleration filters off), -E 10 --domE 10 (every hit at E <= 10;
# no --cut_ga), on results/step36b/chunks/C.<task>.faa (distinct member sequences of the 4,509 recurrent EU families,
# scripts/step36b_prep.py). Outputs results/step36b/hmm/C.<task>.{tblout,domtblout}; refuses to overwrite.
module load hmmer/3.4
module load db-pfam
module load workspace/scratch
W=/bigdata/stajichlab/lshad003/wf_protein_space; D=$W/results/step36b
I=$D/chunks/C.$SLURM_ARRAY_TASK_ID.faa; O=$D/hmm/C.$SLURM_ARRAY_TASK_ID
if [ ! -s $I ]; then echo "ABORT: no input $I"; exit 1; fi
if [ -e $O.tblout ]; then echo "REFUSING: $O.tblout exists"; exit 1; fi
mkdir -p $D/hmm
# node-local job scratch from workspace/scratch (2026-10-24: /scratch/$USER is not writable on compute nodes)
[ -n "$SCRATCH" ] || { echo "ABORT: SCRATCH not set"; exit 1; }
L=$SCRATCH/pfam; mkdir -p $L || { echo "ABORT: cannot create $L"; exit 1; }
for F in $PFAM_DB/Pfam-A.hmm*; do cp $F $L/ || { echo "ABORT: staging $F failed"; exit 1; }; done
hostname; date; echo "PFAM_DB $PFAM_DB; sequences $(grep -c '^>' $I)"
hmmsearch --max -E 10 --domE 10 --cpu $SLURM_CPUS_PER_TASK --tblout $O.tblout.tmp --domtblout $O.domtblout.tmp $L/Pfam-A.hmm $I > /dev/null
rc=$?
if [ $rc -eq 0 ]; then mv $O.tblout.tmp $O.tblout; mv $O.domtblout.tmp $O.domtblout; fi
rm -rf $L
echo "exit $rc; hits $(grep -vc '^#' $O.tblout 2>/dev/null)"; date
exit $rc
