#!/bin/bash
#SBATCH -p epyc -c 32 --mem 240G --time 4-00:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step13f.log
# Step 13f: supported-U reps vs NCBI ClusteredNR 20260128, DIAMOND blastp
# very-sensitive. DB HARD-CODED after step13e dbinfo (470,748,714 seqs).
# Loose record E<=1e-3 k=5; the hit criterion is applied later on fields.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
Q=$WPS/results/supported_U_reps.faa
DB=/srv/projects/db/ncbi/diamond/20260128/nr_cluster_seq.dmnd
OUT=$WPS/results/step13_nrclust20260128_hits.tsv

N=$(grep -c '>' $Q 2>/dev/null)
[ "$N" = "586215" ] || { echo "GUARD FAIL: query count $N need 586215"; exit 1; }
module load diamond
diamond version
SEQ=$(diamond dbinfo --db $DB 2>/dev/null | awk '/Sequences/{print $2}')
[ "$SEQ" = "470748714" ] || { echo "GUARD FAIL: dbinfo $SEQ expect 470748714"; exit 1; }

SCR=/scratch/step13f_$SLURM_JOB_ID
mkdir -p $SCR
AV=$(df -BG --output=avail /scratch 2>/dev/null | tail -1 | tr -dc '0-9')
echo "scratch avail: ${AV}G"
if [ -n "$AV" ] && [ "$AV" -ge 220 ] && cp $DB $SCR/db.dmnd 2>/dev/null; then
  LDB=$SCR/db.dmnd; echo "staged to $SCR"
else
  LDB=$DB; rm -f $SCR/db.dmnd; echo "reading in place, single job"
fi

diamond blastp --query $Q --db $LDB --threads 32 -b 8 -c 1 --very-sensitive \
  -e 1e-3 -k 5 \
  --outfmt 6 qseqid sseqid pident length qlen evalue bitscore qcovhsp scovhsp stitle \
  --tmpdir $SCR --out $OUT 2> $WPS/logs/step13f.diamond.err
echo "queries with any hit at E<=1e-5: $(awk -F'\t' '$6<=1e-5{print $1}' $OUT | sort -u | wc -l) / 586215"
rm -rf $SCR
echo "DONE step13f"
