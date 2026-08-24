#!/bin/bash
#SBATCH -p stajichlab -c 16 --mem 64G --time 24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step16m_pilot.log
# Step 16m: pilot mapping, ONE sample (smallest R1). Counts in one pass under
# two stated criteria: primary mapped reads (no secondary, no supplementary,
# R1+R2 as separate reads) and the MAPQ>=10 subset. No BAM written anywhere.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
REF=$WPS/catalog/db/LsPS_CDS_95_supported.fasta
mkdir -p $WPS/results/counts $WPS/results/mapping_summary

OKJ=$(grep -c 'DONE step16j' $WPS/logs/step16j.log 2>/dev/null)
[ "$OKJ" = "1" ] || { echo "GUARD FAIL: index job not finished clean"; exit 1; }

S=""; R1=""; R2=""; MIN=0
while IFS=$'\t' read -r st r1 r2; do
  [ "$st" = "stem" ] && continue
  SZ=$(stat -c %s "$r1" 2>/dev/null) || continue
  if [ -z "$S" ] || [ "$SZ" -lt "$MIN" ]; then S=$st; R1=$r1; R2=$r2; MIN=$SZ; fi
done < $WPS/results/reads_89.tsv
echo "pilot sample: $S, R1 bytes: $MIN"

module load bwa-mem2 samtools
bwa-mem2 mem -t 16 $REF "$R1" "$R2" 2> $WPS/logs/step16m_pilot.bwa.err \
| samtools view -F 0x900 - \
| awk -v CF="$WPS/results/counts/$S.counts.tsv" -v S="$S" '
    {tot++; if (int($2/4)%2==1) unm++; else {a[$3]++; if ($5>=10) q[$3]++}}
    END{for (g in a) printf "%s\t%d\t%d\n", g, a[g], q[g] > CF;
        printf "%s\t%d\t%d\t%.4f\n", S, tot, tot-unm, (tot-unm)/tot}' \
  > $WPS/results/mapping_summary/$S.summary.tsv
gzip -f $WPS/results/counts/$S.counts.tsv
echo "== stem, records, mapped_primary, fraction =="
cat $WPS/results/mapping_summary/$S.summary.tsv
echo "genes with nonzero count: $(zcat $WPS/results/counts/$S.counts.tsv.gz | wc -l)"
echo "DONE step16m pilot"
