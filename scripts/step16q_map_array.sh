#!/bin/bash
#SBATCH -p stajichlab -c 16 --mem 64G --time 12:00:00
#SBATCH --array=1-89%8
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16q.%a.log
# Step 16q: abundance mapping, all 89, supported-CDS reference. Counts both
# criteria per gene: col2 primary mapped, col3 MAPQ>=10 subset. No BAM.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
REF=$WPS/catalog/db/LsPS_CDS_95_supported.fasta
[ -s $REF.bwt.2bit.64 ] || { echo "GUARD FAIL: index missing"; exit 1; }
mkdir -p $WPS/results/counts $WPS/results/mapping_summary

LINE=$((SLURM_ARRAY_TASK_ID+1))
ROW=$(sed -n "${LINE}p" $WPS/results/reads_89.tsv)
IFS=$'\t' read -r S R1 R2 <<< "$ROW"
[ -n "$S" ] || { echo "no row $LINE"; exit 1; }
if [ -s $WPS/results/counts/$S.counts.tsv.gz ] && [ -s $WPS/results/mapping_summary/$S.summary.tsv ]; then
  echo "SKIP $S, outputs exist"; exit 0
fi
[ -s "$R1" ] && [ -s "$R2" ] || { echo "MISSING READS $S"; exit 1; }

module load bwa-mem2/2.3 samtools
echo "sample $S"
bwa-mem2 mem -t 16 $REF "$R1" "$R2" 2> $WPS/logs/step16q.$S.bwa.err \
| samtools view -F 0x900 - \
| awk -v CF="$WPS/results/counts/$S.counts.tsv" -v S="$S" '
    {tot++; if (int($2/4)%2==1) unm++; else {a[$3]++; if ($5>=10) q[$3]++}}
    END{for (g in a) printf "%s\t%d\t%d\n", g, a[g], q[g] > CF;
        printf "%s\t%d\t%d\t%.4f\n", S, tot, tot-unm, (tot-unm)/tot}' \
  > $WPS/results/mapping_summary/$S.summary.tsv
gzip -f $WPS/results/counts/$S.counts.tsv
cat $WPS/results/mapping_summary/$S.summary.tsv
echo "DONE $S"
