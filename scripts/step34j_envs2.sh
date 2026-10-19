#!/bin/bash
#SBATCH -J step34j_envs
#SBATCH -p epyc
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_envs2.%j.log
# Step 34j: second attempt at the new environments (step34j_envs.sh was killed on the login node: memory limit).
# Same packages; mamba with a private package cache on node scratch. The Merizo clone from the first attempt
# (software/Merizo, commit 41d12fb84e6e8fdb586c2c859d12161dc7bb5bfd) is kept. Refuses to overwrite envs or databases.
C=/bigdata/stajichlab/lshad003/condaenvs
S=/bigdata/stajichlab/lshad003/software
MM=/opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/condabin/mamba
for p in $C/genomad $C/defensefinder $C/merizo $S/genomad_db $S/defensefinder_models; do
  if [ -e $p ]; then echo "REFUSING: $p exists"; exit 1; fi
done
module load workspace/scratch
export CONDA_PKGS_DIRS=$SCRATCH/pkgs; mkdir -p $CONDA_PKGS_DIRS
hostname; date
echo "== genomad"
$MM create -y -q -p $C/genomad -c conda-forge -c bioconda --override-channels genomad
$C/genomad/bin/genomad --version
$C/genomad/bin/genomad download-database $S
cat $S/genomad_db/version.txt; echo
echo "== defensefinder"
$MM create -y -q -p $C/defensefinder -c conda-forge -c bioconda --override-channels python=3.10 hmmer pip
$C/defensefinder/bin/pip install -q mdmparis-defense-finder
$C/defensefinder/bin/defense-finder --version
$C/defensefinder/bin/defense-finder update --models-dir $S/defensefinder_models
ls $S/defensefinder_models
echo "== merizo"
$MM create -y -q -p $C/merizo -c conda-forge --override-channels python=3.10 pip
$C/merizo/bin/pip install -q torch==2.1.2 --index-url https://download.pytorch.org/whl/cu118
$C/merizo/bin/pip install -q numpy==1.26.4 einops scipy networkx rotary_embedding_torch natsort matplotlib
$C/merizo/bin/python -c "import torch; print('torch', torch.__version__, torch.version.cuda)"
date
