#!/bin/bash
# Step 34j: new environments for geNomad, DefenseFinder and Merizo (never rf_py39). Refuses to overwrite.
# geNomad: bioconda genomad, database via `genomad download-database` into software/genomad_db.
# DefenseFinder: pip mdmparis-defense-finder with bioconda hmmer; models via `defense-finder update` into
#   software/defensefinder_models. Merizo: git clone psipred/Merizo into software/Merizo (commit printed), env with
#   pytorch (pip, CUDA 11.8 wheels; also runs on CPU).
C=/bigdata/stajichlab/lshad003/condaenvs
S=/bigdata/stajichlab/lshad003/software
MM=/opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/condabin/mamba
for p in $C/genomad $C/defensefinder $C/merizo $S/genomad_db $S/defensefinder_models $S/Merizo; do
  if [ -e $p ]; then echo "REFUSING: $p exists"; exit 1; fi
done
date
echo "== genomad"
$MM create -y -p $C/genomad -c conda-forge -c bioconda genomad
$C/genomad/bin/genomad --version
$C/genomad/bin/genomad download-database $S && mv $S/genomad_db $S/genomad_db 2>/dev/null; ls $S/genomad_db | head; cat $S/genomad_db/version.txt
echo "== defensefinder"
$MM create -y -p $C/defensefinder -c conda-forge -c bioconda python=3.10 hmmer pip
$C/defensefinder/bin/pip install mdmparis-defense-finder
$C/defensefinder/bin/defense-finder --version
$C/defensefinder/bin/defense-finder update --models-dir $S/defensefinder_models
ls $S/defensefinder_models
echo "== merizo"
git clone https://github.com/psipred/Merizo.git $S/Merizo && git -C $S/Merizo log -1 --format="commit %H %ad"
$MM create -y -p $C/merizo -c conda-forge python=3.10 pip
$C/merizo/bin/pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cu118
$C/merizo/bin/pip install numpy==1.26.4 einops scipy networkx rotary_embedding_torch natsort matplotlib
$C/merizo/bin/python -c "import torch; print('torch', torch.__version__, torch.version.cuda)"
ls $S/Merizo
date
