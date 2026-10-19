#!/bin/bash
#SBATCH -J step34j_envs3
#SBATCH -p epyc
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_envs3.%j.log
# Step 34j: finish the step34j_envs2.sh setup.
# geNomad: bioconda resolved 1.13.0, which is not in the geNomad database compatibility table
#   (portal.nersc.gov/genomad/__data__/releases.txt: db v1.9 for packages 1.10-1.12), so download-database refused.
#   Pin genomad to the newest 1.12 build in the condaenvs/genomad env, then download db v1.9 into software/genomad_db.
# DefenseFinder: `defense-finder update` failed only because hmmsearch was not on PATH (env not activated); rerun
#   with the env bin directory on PATH. Refuses to overwrite the databases.
C=/bigdata/stajichlab/lshad003/condaenvs
S=/bigdata/stajichlab/lshad003/software
MM=/opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/condabin/mamba
if [ -e $S/genomad_db ]; then echo "REFUSING: $S/genomad_db exists"; exit 1; fi
module load workspace/scratch
export CONDA_PKGS_DIRS=$SCRATCH/pkgs; mkdir -p $CONDA_PKGS_DIRS
hostname; date
$MM install -y -q -p $C/genomad -c conda-forge -c bioconda --override-channels "genomad>=1.12,<1.13"
$C/genomad/bin/genomad --version
PATH=$C/genomad/bin:$PATH genomad download-database $S
cat $S/genomad_db/version.txt; echo
ls -A $S/defensefinder_models
PATH=$C/defensefinder/bin:$PATH defense-finder update --models-dir $S/defensefinder_models
ls $S/defensefinder_models; cat $S/defensefinder_models/*/metadata.yml 2>/dev/null | grep -i "vers" | head
$C/defensefinder/bin/pip show mdmparis-defense-finder | grep -i version
date
