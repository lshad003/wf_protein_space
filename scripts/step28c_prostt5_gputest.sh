#!/bin/bash
#SBATCH -p gpu --gres=gpu:a100:1 -c 8 --mem 32G --time 2:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28c.%j.log
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
M=/bigdata/stajichlab/lshad003/software/foldseek_gpu/weights/prostt5
W=/bigdata/stajichlab/lshad003/wf_protein_space
T=$W/tmp_step28c; mkdir -p $T
awk '/^>/{n++} n<=100' $W/results/unk_filtered_aplakidou.faa > $T/test100.faa
echo "test proteins: $(grep -c '^>' $T/test100.faa)"
nvidia-smi --query-gpu=name,memory.total --format=csv
date
$F createdb $T/test100.faa $T/test100_3di --prostt5-model $M --gpu 1
date
ls -la $T
echo "DONE step28c"
