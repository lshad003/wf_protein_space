#!/bin/bash
# Step 16n: 16j verdict after its purge from the queue; index files; queue; pfam.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space

echo "=== 1. step16j final state, job 27713810 ==="
sacct -j 27713810 -n -X --format=State,ExitCode,Elapsed,MaxRSS
tail -6 $BASE/logs/step16j.log 2>/dev/null

echo ""
echo "=== 2. index files on disk ==="
ls -lh $BASE/catalog/db | grep LsPS_CDS_95_supported

echo ""
echo "=== 3. queue: union job 27713842 should sit PENDING with a BeginTime ==="
squeue -u lshad003 -o "%.10i %.9P %.12j %.2t %.10M %.20S %R" | head -15

echo ""
echo "=== 4. pfam, timestamped ==="
date
ls $BASE/results/pfam_hs | grep -c '\.done$'
