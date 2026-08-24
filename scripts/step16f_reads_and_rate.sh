#!/bin/bash
# Step 16f: fastq layout one level down per cohort, Pfam rate from .done mtimes.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal

for F in wf22_stems_44.txt wf23_stems.txt wf24_stems.txt; do
  S=$(head -1 $BASE/metadata/$F)
  echo "=== $S ==="
  ls $FEC/results/$S/qc 2>/dev/null
  D=$(ls -d $FEC/results/$S/qc/*/ 2>/dev/null | head -1)
  if [ -n "$D" ]; then echo "--- inside $D ---"; ls -l $D | head -8; fi
done

echo ""
echo "=== pfam rate from newest 100 .done mtimes ==="
date
DONE=$(ls $BASE/results/pfam_hs | grep -c '\.done$')
echo "done: $DONE / 619"
ls -lt --time-style=+%s $BASE/results/pfam_hs/*.done | head -100 | awk -v d=$DONE '
  {t[NR]=$6}
  END{dt=t[1]-t[NR];
      if(dt>0){r=(NR-1)*3600/dt;
        printf "newest %d span %.2f h, rate %.1f chunks per h, eta for %d left: %.1f h\n",
               NR, dt/3600, r, 619-d, (619-d)/r}}'
