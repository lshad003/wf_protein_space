#!/bin/bash
#SBATCH -p epyc -c 32 --mem 128G --time 2-00:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step13d.log
# Step 13d: supported-U reps vs UniRef90, DIAMOND very-sensitive. Loose record
# (E<=1e-3, k=5); the hit criterion is applied later on recorded fields.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  $WPS/scripts/step13d_extract_u.py || exit 1

NEW=$(ls /srv/projects/db/Uniprot/2025_03/*.dmnd 2>/dev/null | head -1)
DB=${NEW:-/srv/projects/db/Uniprot/uniref90.dmnd}
echo "database: $DB"
module load diamond 2>/dev/null
command -v diamond > /dev/null 2>&1 || { echo "no diamond module:"; module -t avail 2>&1 | grep -i diamond; exit 1; }
diamond version

SCR=/scratch/step13d_$SLURM_JOB_ID
mkdir -p $SCR
if cp $DB $SCR/db.dmnd 2>/dev/null; then LDB=$SCR/db.dmnd; echo "staged to $SCR"
else LDB=$DB; echo "scratch copy failed, reading in place"; fi
diamond dbinfo --db $LDB > /dev/null 2>&1 || { echo "dmnd unreadable by this diamond"; rm -rf $SCR; exit 1; }

diamond blastp --query $WPS/results/supported_U_reps.faa --db $LDB \
  --threads 32 -b 8 -c 1 --very-sensitive -e 1e-3 -k 5 \
  --outfmt 6 qseqid sseqid pident length qlen evalue bitscore qcovhsp scovhsp stitle \
  --out $WPS/results/step13_uniref90_hits.tsv 2> $WPS/logs/step13d.diamond.err
echo "queries with any hit at E<=1e-5: $(awk -F'\t' '$6<=1e-5{print $1}' $WPS/results/step13_uniref90_hits.tsv | sort -u | wc -l) / 586215"
rm -rf $SCR
echo "DONE step13d"
