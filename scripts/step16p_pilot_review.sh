#!/bin/bash
# Step 16p: pilot criteria comparison, symlink confirmation, MaxRSS.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
S=UHM17.35739

echo "=== 1. pilot memory and steps ==="
sacct -j 27713872 -n --format=JobID%20,MaxRSS,Elapsed | head -4

echo ""
echo "=== 2. the 85-byte R1: symlink check ==="
R1=$(awk -F'\t' -v s=$S '$1==s{print $2}' $BASE/results/reads_89.tsv)
ls -l "$R1"

echo ""
echo "=== 3. criteria comparison from the pilot counts ==="
zcat $BASE/results/counts/$S.counts.tsv.gz | awk -F'\t' '
  {p+=$2; q+=$3; if($2>0)np++; if($3>0)nq++}
  END{printf "primary sum      %d  (expect 134467558)\n", p;
      printf "mapq10 sum       %d  (%.2f%% of primary)\n", q, 100*q/p;
      printf "genes nonzero    primary %d, mapq10 %d\n", np, nq}'
