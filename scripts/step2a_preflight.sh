#!/bin/bash
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
echo "== mmseqs2 modules available =="
module avail mmseqs2 2>&1 | head -10
echo
echo "== version that loads by default =="
module load mmseqs2 && mmseqs version
echo
echo "== v1 database dbtype / version markers =="
ls -l $V1/db/LsFMGC_AA.dbtype $V1/db/LsFMGC_AA.index $V1/db/LsFMGC_AA.lookup
head -c 20 $V1/db/LsFMGC_AA.dbtype | od -c | head -2
ls $V1/db | grep -i "version\|\.sh" | head
echo
echo "== verify the 5,055,108 claim =="
grep -c "^>" $V1/db/LsFMGC_AA_95_rep.fasta
echo
echo "== rep counts for the other tiers =="
for T in 100 90 50; do
  echo -n "tier $T: "; grep -c "^>" $V1/db/LsFMGC_AA_${T}_rep.fasta
done
echo
echo "== disk space in my space =="
df -h /bigdata/stajichlab/lshad003 | tail -2
du -sh /bigdata/stajichlab/lshad003/wf_protein_space/catalog 2>/dev/null
