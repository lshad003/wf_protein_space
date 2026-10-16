#!/bin/bash
#SBATCH -J step33g_control
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 6:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33g_control.%j.log
# Step 33g control: minimap2 2.30 -x asm10 index of D. melanogaster dm6 (catalog/fly/, see SOURCE.txt) ->
# catalog/fly/dm6.mmi; 5,000 simulated fly fragments with 1% substitutions (scripts/step33g_sim.py) mapped with
# -x asm10 --secondary=no -c; scripts/step33g_sim_summary.py. Exits with the summary's code (4 = fewer than 90% of
# fragments recovered at >= 95%; the step33g chain stops).
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
F=$W/catalog/fly; D=$W/results/step33g_sim
for p in $F/dm6.mmi $D/sim5000.fly.paf; do if [ -e $p ]; then echo "REFUSING: $p exists"; exit 1; fi; done
hostname; date
module load minimap2/2.30
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 -d $F/dm6.mmi $F/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz || { echo "exit 1 (index)"; exit 1; }
$PY $W/scripts/step33g_sim.py || { echo "exit 1 (sim)"; exit 1; }
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no -c $F/dm6.mmi $D/sim5000.fa > $D/sim5000.fly.paf || { echo "exit 1 (map)"; exit 1; }
$PY $W/scripts/step33g_sim_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
