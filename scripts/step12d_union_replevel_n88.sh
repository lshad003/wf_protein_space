#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 16G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step12d_n88.log
# Step 12d n88: rep-level eggNOG+Pfam union on the 88-metagenome catalog.
# Annotation is per representative and unchanged; hit lists from the 89 run
# (results/eggnog_hit_ids.txt, pfam_hit_ids.txt, re-derived here and checked)
# are restricted to results/reps_95_n88.txt. Supported set is
# results/supported_reps_95_n88.txt. Criteria as step12d: eggNOG hit = any
# emapper assignment; Pfam hit = tblout presence under --cut_ga.
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
PF=$WPS/results/pfam_hs
R=$WPS/results
export LC_ALL=C

DONE=$(ls $PF | grep -c '\.done$')
[ "$DONE" = "619" ] || { echo "GUARD FAIL: done=$DONE need 619"; exit 1; }

cat $PF/*.tblout | grep -v '^#' | awk '{print $1}' | sort -u -S 8G -T $WPS/tmp > $R/pfam_hit_ids_n89check.tmp
cmp -s $R/pfam_hit_ids_n89check.tmp $R/pfam_hit_ids.txt || { echo "CROSS-CHECK FAIL: pfam ids differ from pfam_hit_ids.txt"; exit 1; }
rm $R/pfam_hit_ids_n89check.tmp
NE89=$(wc -l < $R/eggnog_hit_ids.txt)
[ "$NE89" = "2984013" ] || { echo "CROSS-CHECK FAIL: eggnog_hit_ids.txt has $NE89"; exit 1; }
echo "89 inputs verified: pfam_hit_ids.txt reproduced, eggnog_hit_ids.txt 2984013"

sort -c $R/reps_95_n88.txt && sort -c $R/supported_reps_95_n88.txt || { echo "inputs not LC_ALL=C sorted"; exit 1; }
comm -12 $R/eggnog_hit_ids.txt $R/reps_95_n88.txt > $R/eggnog_hit_ids_n88.txt
comm -12 $R/pfam_hit_ids.txt   $R/reps_95_n88.txt > $R/pfam_hit_ids_n88.txt
sort -u -S 8G -T $WPS/tmp $R/eggnog_hit_ids_n88.txt $R/pfam_hit_ids_n88.txt > $R/known_union_ids_n88.txt
NE=$(wc -l < $R/eggnog_hit_ids_n88.txt); NP=$(wc -l < $R/pfam_hit_ids_n88.txt)
UNION=$(wc -l < $R/known_union_ids_n88.txt)
BOTH=$(comm -12 $R/eggnog_hit_ids_n88.txt $R/pfam_hit_ids_n88.txt | wc -l)
SE=$(comm -12 $R/eggnog_hit_ids_n88.txt $R/supported_reps_95_n88.txt | wc -l)
SP=$(comm -12 $R/pfam_hit_ids_n88.txt  $R/supported_reps_95_n88.txt | wc -l)
SB=$(comm -12 $R/eggnog_hit_ids_n88.txt $R/pfam_hit_ids_n88.txt | comm -12 - $R/supported_reps_95_n88.txt | wc -l)
SU=$(comm -12 $R/known_union_ids_n88.txt $R/supported_reps_95_n88.txt | wc -l)
T=$(wc -l < $R/reps_95_n88.txt); TS=$(wc -l < $R/supported_reps_95_n88.txt)
NOTIN=$(comm -23 $R/supported_reps_95_n88.txt $R/reps_95_n88.txt | wc -l)
echo "supported reps not in reps_95_n88: $NOTIN (must be 0)"; [ "$NOTIN" = "0" ] || exit 1
{
echo "# n88. eggNOG hit = any emapper assignment; Pfam hit = tblout under --cut_ga"
printf "scope\ttotal\teggnog\tpfam\tboth\tunion\tneither\n"
printf "all_reps\t%d\t%d\t%d\t%d\t%d\t%d\n" $T $NE $NP $BOTH $UNION $((T-UNION))
printf "supported\t%d\t%d\t%d\t%d\t%d\t%d\n" $TS $SE $SP $SB $SU $((TS-SU))
} > $R/union_rep_level_n88.tsv
cat $R/union_rep_level_n88.tsv
awk -F'\t' '$1!~/^#/ && $1!="scope"{printf "%s eggnog %.2f%% pfam %.2f%% union %.2f%% neither %.2f%%\n",$1,100*$3/$2,100*$4/$2,100*$6/$2,100*$7/$2}' $R/union_rep_level_n88.tsv
echo "DONE step12d_n88"
