#!/bin/bash
#SBATCH -p epyc -c 2 --mem 8G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28a2.%j.log
D=/bigdata/stajichlab/lshad003/software/foldseek_gpu
cd $D
$D/foldseek/bin/foldseek version
for i in 1 2 3 4 5; do
  $D/foldseek/bin/foldseek databases ProstT5 $D/weights/prostt5 $D/tmp && break
  sleep 60
done
ls -la $D/weights
echo "DONE step28a2"
