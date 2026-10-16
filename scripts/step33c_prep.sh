#!/bin/bash
#SBATCH -J step33c_prep
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 64G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33c_prep.%j.log
# Step 33c prep: contig lists (scripts/step33c_prep.py), then a minimap2 2.30 -x asm10 index of the six Basidiobolus
# genomes (catalog/basid/basid6.mmi). Exits non-zero on any failure.
W=/bigdata/stajichlab/lshad003/wf_protein_space
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step33c_prep.py || { echo "exit 1 (prep.py)"; exit 1; }
module load minimap2/2.30
if [ -e $W/catalog/basid/basid6.mmi ]; then echo "REFUSING: basid6.mmi exists"; exit 1; fi
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 -d $W/catalog/basid/basid6.mmi $W/catalog/basid/basid6.fa
rc=$?
echo "exit $rc"; date
exit $rc
