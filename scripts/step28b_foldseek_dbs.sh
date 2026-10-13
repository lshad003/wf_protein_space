#!/bin/bash
#SBATCH -p epyc -c 4 --mem 16G --time 24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28b.%j.log
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
D=/bigdata/stajichlab/lshad003/software/foldseek_dbs
mkdir -p $D/tmp && cd $D
$F version
for db in PDB CATH50; do
  for i in 1 2 3 4 5; do
    $F databases $db $D/$db $D/tmp && break
    sleep 120
  done
done
ls -la $D
du -sh $D
echo "DONE step28b"
