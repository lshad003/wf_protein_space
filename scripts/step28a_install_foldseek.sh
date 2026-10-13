#!/bin/bash
#SBATCH -p stajichlab -c 2 --mem 8G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28a.%j.log
D=/bigdata/stajichlab/lshad003/software/foldseek_gpu
mkdir -p $D/weights $D/tmp && cd $D
for i in 1 2 3 4 5; do
  curl -L -C - -o foldseek-linux-gpu.tar.gz https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz && break
  sleep 60
done
tar xzf foldseek-linux-gpu.tar.gz
$D/foldseek/bin/foldseek version
for i in 1 2 3 4 5; do
  $D/foldseek/bin/foldseek databases ProstT5 $D/weights/prostt5 $D/tmp && break
  sleep 60
done
ls -la $D/weights
df -h /bigdata/stajichlab/lshad003 | tail -1
echo "DONE step28a"
