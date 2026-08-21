#!/bin/bash
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
FEC=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal
echo "== v1/scripts =="
ls "$V1/scripts"
echo "== v1/gene_prediction (was empty-looking before, confirm) =="
ls -a "$V1/gene_prediction"
echo "== Fecal top level =="
ls "$FEC"
echo "== Fecal/Proteins if present (old INDIR comment in 00_setup) =="
ls "$FEC/Proteins" 2>/dev/null || echo "(no Proteins dir)"
echo "== Fecal scripts-like dirs =="
ls "$FEC/scripts" 2>/dev/null || echo "(no scripts dir)"
ls "$FEC/pipeline" 2>/dev/null || echo "(no pipeline dir)"
echo "== prodigal version recorded in an existing gff =="
zcat "$V1/by_assembly/UHM574.41000.gff.gz" 2>/dev/null | head -3
echo "== slurm log header from the original prediction run =="
head -20 "$V1/by_assembly/slurm-15645413.out"
echo "== contigs exist for the 14 missing stems? =="
for s in UHM585.41009 UHM617.41023 UHM624.41024 UHM626.41025 UHM628.41026 UHM632.41027 UHM633.41028 UHM634.41029 UHM636.41030 UHM644.41031 UHM648.41909 UHM649.41032 UHM656.41033 UHM666.41034; do
  ls "$FEC/results/$s/megahit" 2>/dev/null | head -2 | sed "s|^|$s: |" || echo "$s: NO megahit dir"
done
