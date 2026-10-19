#!/bin/bash
#SBATCH -J step34j_fs
#SBATCH -p epyc
#SBATCH -c 64
#SBATCH --mem 480G
#SBATCH -t 7-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_b3_foldseek.%x.%j.log
# Step 34j B3: Foldseek (GPU build 463739e0..., run on CPU as Step 29e) of the ESMFold models (results/step34j/
# esmfold_pdb; 4,250 EU + 4,509 controls) and the Merizo domains of models >= 200 aa (results/step34j/domains)
# against one database, given as $1: PDB, CATH50, AFDB50 or ESMAtlas30 (software/foldseek_dbs).
# Settings: Step 29e, -e 10 and the Step 29e output fields; hits are filtered at E <= 1e-3 downstream.
#   PDB, CATH50: --exhaustive-search 1 (exactly Step 29e).
#   AFDB50, ESMAtlas30: prefilter mode (-s 9.5 --max-seqs 1000); exhaustive alignment against 53M/37M entries for
#   ~10,000 queries was estimated at ~600 node-hours. This is the only deviation from Step 29e.
# The database is staged to node-local scratch. Output: results/step34j/foldseek/<db>.m8. Refuses to overwrite.
DBN=$1
W=/bigdata/stajichlab/lshad003/wf_protein_space
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
O=$W/results/step34j/foldseek
case $DBN in
  PDB|CATH50) MODE="--exhaustive-search 1" ;;
  AFDB50|ESMAtlas30) MODE="-s 9.5 --max-seqs 1000" ;;
  *) echo "ABORT: unknown database $DBN"; exit 1 ;;
esac
if [ -e $O/$DBN.m8 ]; then echo "REFUSING: $O/$DBN.m8 exists"; exit 1; fi
if [ ! -e $DB/$DBN.dbtype ]; then echo "ABORT: $DB/$DBN missing"; exit 1; fi
mkdir -p $O
hostname; date; $F version; echo "db $DBN mode $MODE"
module load workspace/scratch
T=$SCRATCH/step34j_fs_$DBN; mkdir -p $T/db
for f in $DB/$DBN $DB/${DBN}_*  $DB/$DBN.*; do [ -e "$f" ] && cp $f $T/db/; done
du -sh $T/db; date
# 2026-10-19 fix: createdb accepts one directory or one list file, not two directories (jobs 30050462-30050465 failed
# with "Only one directory or tsv file ... can be given"); queries are now passed as a list of model and domain PDBs.
# Guards: domains must be non-empty, UHM125.35769__k141_61975_1 (results/step34j/esmfold_excluded.tsv) must be absent,
# createdb must load every listed file, and search must exit 0 before convertalis.
nm=$(ls $W/results/step34j/esmfold_pdb | grep -c '\.pdb$'); nd=$(ls $W/results/step34j/domains | grep -c '\.pdb$')
if [ $nd -eq 0 ]; then echo "ABORT: results/step34j/domains has no PDBs"; rm -rf $T; exit 1; fi
ls $W/results/step34j/esmfold_pdb/*.pdb > $T/q_files.tsv
ls $W/results/step34j/domains/*.pdb >> $T/q_files.tsv
if grep -q "UHM125.35769__k141_61975_1" $T/q_files.tsv; then echo "ABORT: excluded UHM125.35769__k141_61975_1 in query list"; rm -rf $T; exit 1; fi
echo "query files: models $nm, domains $nd, listed $(wc -l < $T/q_files.tsv)"
$F createdb $T/q_files.tsv $T/q --threads $SLURM_CPUS_PER_TASK
if [ ! -e $T/q.index ]; then echo "ABORT: createdb wrote no query database"; rm -rf $T; exit 1; fi
echo "query entries $(wc -l < $T/q.index)"
if [ $(wc -l < $T/q.index) -ne $((nm + nd)) ]; then echo "ABORT: query entries differ from $((nm + nd)) files"; rm -rf $T; exit 1; fi
$F search $T/q $T/db/$DBN $T/res $T/tmp -e 10 $MODE --threads $SLURM_CPUS_PER_TASK
rc=$?
if [ $rc -ne 0 ]; then echo "ABORT: search exit $rc"; rm -rf $T; exit 1; fi
$F convertalis $T/q $T/db/$DBN $T/res $T/out.m8 --threads $SLURM_CPUS_PER_TASK \
  --format-output query,target,fident,alnlen,qlen,tlen,evalue,bits,prob,alntmscore,qtmscore,ttmscore,lddt,qcov,tcov
rc=$?
echo "search exit 0 convertalis exit $rc"
if [ $rc -ne 0 ] || [ ! -e $T/out.m8 ]; then echo "ABORT: convertalis exit $rc or no output"; rm -rf $T; exit 1; fi
cp $T/out.m8 $O/$DBN.m8
rc=$?
if [ $rc -ne 0 ]; then echo "ABORT: copy to $O/$DBN.m8 failed"; rm -rf $T; exit 1; fi
echo "$DBN rows $(wc -l < $O/$DBN.m8); queries with a hit $(cut -f1 $O/$DBN.m8 | sort -u | wc -l); with E <= 1e-3 $(awk -F'\t' '$7<=1e-3{print $1}' $O/$DBN.m8 | sort -u | wc -l)"
rm -rf $T; date
# 2026-10-19: exit with the main commands' status so afterok dependencies see failures (CLAUDE.md rule)
exit $rc
