#!/bin/bash
#SBATCH -J step34j_mge
#SBATCH -p epyc
#SBATCH -c 32
#SBATCH --mem 128G
#SBATCH -t 7-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_a3_mge.%j.log
# Step 34j A3: geNomad (condaenvs/genomad, database software/genomad_db staged to node scratch; end-to-end, default
# settings) on results/step34j/mge/contigs.fna, and DefenseFinder (condaenvs/defensefinder, models
# software/defensefinder_models; --db-type gembase, default coverage and cut_ga) on results/step34j/mge/proteins.faa.
# Inputs from scripts/step34j_a3_prep.py (run here). Outputs: results/step34j/mge/genomad/, results/step34j/mge/defensefinder/.
W=/bigdata/stajichlab/lshad003/wf_protein_space
C=/bigdata/stajichlab/lshad003/condaenvs
S=/bigdata/stajichlab/lshad003/software
O=$W/results/step34j/mge
hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step34j_a3_prep.py || { echo "ABORT: prep"; exit 1; }
module load workspace/scratch
T=$SCRATCH/step34j_mge; mkdir -p $T
cp -r $S/genomad_db $T/
echo "== geNomad $($C/genomad/bin/genomad --version) db $(cat $T/genomad_db/version.txt)"; date
PATH=$C/genomad/bin:$PATH genomad end-to-end --cleanup --threads $SLURM_CPUS_PER_TASK $O/contigs.fna $O/genomad $T/genomad_db
echo "genomad exit $?"; date
echo "== DefenseFinder $($C/defensefinder/bin/pip show mdmparis-defense-finder | grep -i ^version)"
PATH=$C/defensefinder/bin:$PATH defense-finder run --db-type gembase --models-dir $S/defensefinder_models \
  -w $SLURM_CPUS_PER_TASK -o $O/defensefinder $O/proteins.faa
echo "defensefinder exit $?"; date
ls $O/genomad $O/defensefinder
rm -rf $T
