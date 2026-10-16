#!/bin/bash
#SBATCH -J step33d_poscontrol
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem 64G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step33d_poscontrol.%j.log
# Step 33d: host-detection positive control. scripts/step33d_sim.py (5,000 simulated frog fragments, 1% substitutions),
# minimap2 2.30 -x asm10 --secondary=no against catalog/host/aRanSyl1.mmi (staged to node-local $SCRATCH, as step33c),
# then scripts/step33d_summary.py. Exits non-zero on any failure.
W=/bigdata/stajichlab/lshad003/wf_protein_space
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
D=$W/results/step33d_sim
if [ -e $D/sim5000.host.paf ]; then echo "REFUSING: $D/sim5000.host.paf exists"; exit 1; fi
hostname; date
$PY $W/scripts/step33d_sim.py || { echo "exit 1 (sim)"; exit 1; }
module load minimap2/2.30
module load workspace/scratch
cp $W/catalog/host/aRanSyl1.mmi $SCRATCH/host.mmi || { echo "ABORT: staging failed"; exit 1; }
minimap2 -t $SLURM_CPUS_PER_TASK -x asm10 --secondary=no $SCRATCH/host.mmi $D/sim5000.fa > $SCRATCH/sim.paf
rc=$?
echo "minimap2 exit $rc"
if [ $rc -ne 0 ]; then exit 1; fi
cp $SCRATCH/sim.paf $D/sim5000.host.paf
$PY $W/scripts/step33d_summary.py
rc=$?
echo "exit $rc"; date
exit $rc
