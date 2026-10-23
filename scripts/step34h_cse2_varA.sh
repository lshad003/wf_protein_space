#!/bin/bash
#SBATCH -J step34h_varA
#SBATCH -p epyc -c 32 --mem 128G -t 1-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34h_cse2_varA.%j.log
# Step 34h Part E: ESMFold of the Cse2 variant A rep UHM45.35758__k141_100597_236 with the step29b setup
# (esm 2.0.0 esm-fold, esmfold_3B_v1, TORCH_HOME=/srv/projects/db/esm, --chunk-size 64 --cpu-only, terminal "*"
# removed); model copied to results/esmfold_pdb/UHM45.35758__k141_100597_236.pdb (refuses to overwrite). Foldseek with
# the step29e settings (GPU build commit 463739e run on CPU; PDB and CATH50 from software/foldseek_dbs;
# easy-search --exhaustive-search 1 -e 10; step29e output fields) to results/step34h/foldseek_{PDB,CATH50}.m8.
# Then scripts/fig_render_cse2_varA.py with pymol -cq (8 threads) and scripts/step34h_cse2_varA.py (summary).
W=/bigdata/stajichlab/lshad003/wf_protein_space; D=$W/results/step34h/cse2_varA; R=UHM45.35758__k141_100597_236
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek; DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
for p in $W/results/esmfold_pdb/$R.pdb $W/results/fig_render_cse2_varA.png $D; do [ -e $p ] && { echo "REFUSING: $p exists"; exit 1; }; done
mkdir -p $D; hostname; date
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - $D/query.faa <<'PYEOF'
import sys
want = "UHM45.35758__k141_100597_236"; s = None; h = None; buf = []
for l in open("/bigdata/stajichlab/lshad003/wf_protein_space/results/step34f_cse2/prodigal.faa"):
    if l[0] == ">":
        if h == want: s = "".join(buf)
        h = l[1:].split()[0]; buf = []
    else: buf.append(l.strip())
if h == want: s = "".join(buf)
s = s.rstrip("*"); assert "*" not in s and len(s) == 161
open(sys.argv[1], "w").write(">%s\n%s\n" % (want, s)); print("query length", len(s))
PYEOF
export TORCH_HOME=/srv/projects/db/esm OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
/opt/linux/rocky/8.x/x86_64/pkgs/esm/2.0.0/env/bin/esm-fold -i $D/query.faa -o $D/esmfold --chunk-size 64 --cpu-only 2>&1 | tee $D/esmfold.log
rc=${PIPESTATUS[0]}
if [ $rc -ne 0 ]; then echo "ABORT: esm-fold exit $rc"; exit $rc; fi
[ -s $D/esmfold/$R.pdb ] || { echo "ABORT: no model"; exit 1; }
cp $D/esmfold/$R.pdb $W/results/esmfold_pdb/$R.pdb
module load workspace/scratch; mkdir -p $D/q; cp $D/esmfold/$R.pdb $D/q/
$F version
for db in PDB CATH50; do
  $F easy-search $D/q $DB/$db $D/foldseek_$db.m8 $SCRATCH/tmp_$db --exhaustive-search 1 -e 10 --threads 8 \
    --format-output query,target,fident,alnlen,qlen,tlen,evalue,bits,prob,alntmscore,qtmscore,ttmscore,lddt,qcov,tcov > /dev/null
  echo "$db rows $(wc -l < $D/foldseek_$db.m8)"
done
cd $W/results && OMP_NUM_THREADS=8 /bigdata/stajichlab/lshad003/condaenvs/pymol/bin/pymol -cq $W/scripts/fig_render_cse2_varA.py 2>&1 | tee $D/pymol.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $W/scripts/step34h_cse2_varA.py
rc=$?
echo "exit $rc"; date
# 2026-10-23: esm-fold status taken from PIPESTATUS (tee masked it); exit with the last step's status (CLAUDE.md rule)
exit $rc
