#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 16G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step12d.log
# Step 12d v2: rep-level eggNOG+Pfam union. Replaces v1, whose supported row
# misaligned with its header; v1 never ran past its guard and wrote nothing.
# Criteria: eggNOG hit = any emapper assignment; Pfam hit = tblout presence
# under --cut_ga. Denominators: 6,182,117 reps; 2,069,453 supported.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
PF=$WPS/results/pfam_hs
R=$WPS/results
export LC_ALL=C

DONE=$(ls $PF | grep -c '\.done$')
RUN=$(squeue -u lshad003 -h -n pfamS | wc -l)
if [ "$DONE" != "619" ] || [ "$RUN" != "0" ]; then
  echo "GUARD FAIL: done=$DONE need 619, pfamS in queue=$RUN need 0"; exit 1
fi

cat $PF/*.tblout | grep -v '^#' | awk '{print $1}' | sort -u -S 8G -T $WPS/tmp > $R/pfam_hit_ids.txt
NP=$(wc -l < $R/pfam_hit_ids.txt)
echo "pfam hit reps: $NP"

cat $R/eggnog2/*.emapper.annotations | grep -v '^#' | awk -F'\t' '{print $1}' | sort -u -S 8G -T $WPS/tmp > $R/eggnog_hit_ids.txt
NE=$(wc -l < $R/eggnog_hit_ids.txt)
echo "eggnog hit reps: $NE, expect 2984013"
[ "$NE" = "2984013" ] || { echo "CROSS-CHECK FAIL, stop"; exit 1; }

sort -u -S 8G -T $WPS/tmp $R/eggnog_hit_ids.txt $R/pfam_hit_ids.txt > $R/known_union_ids.txt
UNION=$(wc -l < $R/known_union_ids.txt)
BOTH=$(comm -12 $R/eggnog_hit_ids.txt $R/pfam_hit_ids.txt | wc -l)
SE=$(comm -12 $R/eggnog_hit_ids.txt $R/supported_reps_95.txt | wc -l)
SP=$(comm -12 $R/pfam_hit_ids.txt  $R/supported_reps_95.txt | wc -l)
SB=$(comm -12 $R/eggnog_hit_ids.txt $R/pfam_hit_ids.txt | comm -12 - $R/supported_reps_95.txt | wc -l)
SU=$(comm -12 $R/known_union_ids.txt $R/supported_reps_95.txt | wc -l)

T=6182117; TS=2069453
{
echo "# eggNOG hit = any emapper assignment; Pfam hit = tblout under --cut_ga"
printf "scope\ttotal\teggnog\tpfam\tboth\tunion\tneither\n"
printf "all_reps\t%d\t%d\t%d\t%d\t%d\t%d\n" $T $NE $NP $BOTH $UNION $((T-UNION))
printf "supported\t%d\t%d\t%d\t%d\t%d\t%d\n" $TS $SE $SP $SB $SU $((TS-SU))
} > $R/union_rep_level.tsv
cat $R/union_rep_level.tsv
echo "DONE step12d"
