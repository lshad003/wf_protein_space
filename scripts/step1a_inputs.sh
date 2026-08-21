#!/bin/bash
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
echo "=== 10_predict_by_assembl.sh ==="
cat "$FEC/pipeline/10_predict_by_assembl.sh"
echo
echo "=== megahit dir of a MISSING stem (UHM585.41009) ==="
ls -l "$FEC/results/UHM585.41009/megahit"
echo "--- one level deeper ---"
ls -l "$FEC/results/UHM585.41009/megahit/UHM585.41009_R" 2>/dev/null
echo
echo "=== same for a stem that ALREADY has predictions (UHM574.41000) ==="
ls -l "$FEC/results/UHM574.41000/megahit"
echo
echo "=== scaffolds dir, for comparison ==="
ls -l "$FEC/results/UHM585.41009/scaffolds"
echo
echo "=== prodigal availability ==="
which prodigal || module avail prodigal 2>&1 | head -10
