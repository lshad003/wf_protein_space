#!/bin/bash
#SBATCH -J step33d2_poscontrol
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 64G
#SBATCH -t 6:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33d2_poscontrol.%j.log
# Step 33d v2: the step33d simulated fragments (results/step33d_sim/sim5000.fa, unchanged) mapped with the step33c2
# settings: minimap2 2.30 -x asm10 --secondary=no -c, catalog/host/aRanSyl1.mmi staged to node-local $SCRATCH.
# Then scripts/step33d2_summary.py; exits with its code (4 = fewer than 90% host at >= 95%, positive control failed).
W=/bigdata/stajichlab/lshad003/wf_protein_space
D=$W/results/step33d_sim_v2
if [ -e $D/sim5000.host.paf ]; then echo "REFUSING: $D/sim5000.host.paf exists"; exit 1; fi
mkdir -p $D
hostname; date
module load minimap2/2.30
module load workspace/scratch
cp $W/catalog/host/aRanSyl1.mmi $SCRATCH/host.mmi || { echo "ABORT: staging failed"; exit 1; }
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no -c $SCRATCH/host.mmi $W/results/step33d_sim/sim5000.fa > $SCRATCH/sim.paf
rc=$?
echo "minimap2 exit $rc"
if [ $rc -ne 0 ]; then exit 1; fi
cp $SCRATCH/sim.paf $D/sim5000.host.paf
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step33d2_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
