#!/bin/bash
# Read-only check of results/step34h_member_classes.tsv (one row per strict_final EU 95% cluster).
# Prints: cluster counts under the any and majority rules (Fig. 2e); recurrent families by the highest
# class_any across their clusters (member-validated = EU); member-validated families in which no cluster
# has a distinct eligible sequence beyond its representative (strict_unique_variants <= 1).
W=/bigdata/stajichlab/lshad003/wf_protein_space
awk -F'\t' '
BEGIN{r["EU"]=0; r["GU"]=1; r["KWP"]=2; r["K"]=3}
NR==1{ if($2!="family_rep" || $3!="recurrent_family" || $5!="strict_unique_variants" || $12!="class_any" || $13!="class_majority"){print "ABORT: unexpected header"; bad=1; exit} next }
{a[$12]++; m[$13]++}
$3==1{f=$2; if(!(f in best) || r[$12]>r[best[f]]) best[f]=$12; if($5>1) extra[f]=1}
END{
  if(bad) exit 1
  for(k in a) print "clusters, any rule", k, a[k]
  for(k in m) print "clusters, majority rule", k, m[k]
  for(f in best){n++; c[best[f]]++; if(best[f]=="EU" && !(f in extra)) solo++}
  print "recurrent families", n
  for(k in c) print "recurrent families, highest member class", k, c[k]
  print "member-validated families with no extra member sequence searched", solo+0
}' $W/results/step34h_member_classes.tsv | tee $W/results/check_member_validated.txt
