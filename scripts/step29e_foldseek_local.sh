#!/bin/bash
#SBATCH -J step29e_foldseek
#SBATCH -p epyc
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 4:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step29e_foldseek_local.%j.log
# Step 29e: local Foldseek of the 10 folded no-hit ESMFold models (results/esmfold_peek_v2.tsv: class "folded",
# group NOHIT or NOHIT_LONG), with TM-scores, replacing the by-hand web-server checks (no output files).
# Binary: Foldseek GPU build used in Step 28 (commit 463739e0014a1549a527de589102cde98f802f37), run on CPU.
# Databases: PDB and CATH50 from /bigdata/stajichlab/lshad003/software/foldseek_dbs (Step 28; these carry _ca
# coordinates, so TM-score and lddt are available). No AlphaFoldDB or MGnify Foldseek database exists on the cluster.
# Search: easy-search --exhaustive-search 1 -e 10 (no prefilter, report everything; filter in the summary).
# Summary: top hit by bits per query and database; pass = qTM >= 0.5 and prob >= 0.5.
# Outputs: results/foldseek_folded10/ (query PDBs), results/step29e_foldseek_{PDB,CATH50}.m8,
# results/step29e_foldseek_summary.tsv. Refuses to overwrite; results/foldseek_web4/ is not touched.
W=/bigdata/stajichlab/lshad003/wf_protein_space
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
Q=$W/results/foldseek_folded10
SUM=$W/results/step29e_foldseek_summary.tsv
for p in $Q $SUM $W/results/step29e_foldseek_PDB.m8 $W/results/step29e_foldseek_CATH50.m8; do
  if [ -e $p ]; then echo "REFUSING: $p exists"; exit 1; fi
done
hostname; date; $F version
awk -F'\t' 'NR>1 && $10=="folded" && ($2=="NOHIT" || $2=="NOHIT_LONG") {print $1}' $W/results/esmfold_peek_v2.tsv > $W/results/step29e_query_ids.tmp
n=$(wc -l < $W/results/step29e_query_ids.tmp)
if [ "$n" -ne 10 ]; then echo "ABORT: expected 10 folded no-hit models, got $n"; rm -f $W/results/step29e_query_ids.tmp; exit 1; fi
mkdir -p $Q
while read r; do
  if [ ! -s $W/results/esmfold_pdb/$r.pdb ]; then echo "ABORT: missing $r.pdb"; exit 1; fi
  cp $W/results/esmfold_pdb/$r.pdb $Q/
done < $W/results/step29e_query_ids.tmp
rm -f $W/results/step29e_query_ids.tmp
echo "queries: $(ls $Q | grep -c '\.pdb$')"
module load workspace/scratch
T=$SCRATCH/step29e
mkdir -p $T
for db in PDB CATH50; do
  echo "== $db"; date
  $F easy-search $Q $DB/$db $T/$db.m8 $T/tmp_$db --exhaustive-search 1 -e 10 --threads $SLURM_CPUS_PER_TASK \
    --format-output query,target,fident,alnlen,qlen,tlen,evalue,bits,prob,alntmscore,qtmscore,ttmscore,lddt,qcov,tcov
  echo "exit $?"
  if [ ! -e $T/$db.m8 ]; then echo "ABORT: no output for $db"; rm -rf $T; exit 1; fi
  cp $T/$db.m8 $W/results/step29e_foldseek_$db.m8
  echo "$db: rows $(wc -l < $W/results/step29e_foldseek_$db.m8), queries with a hit $(cut -f1 $W/results/step29e_foldseek_$db.m8 | sort -u | wc -l)"
done
rm -rf $T
$PY - <<'PYEOF'
import os
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
Q = sorted(f[:-4] for f in os.listdir(W + "/results/foldseek_folded10") if f.endswith(".pdb"))
lab = {}
for l in list(open(W + "/results/esmfold_peek_v2.tsv"))[1:]:
    f = l.rstrip("\n").split("\t"); lab[f[0]] = f
def qid(x):  # Foldseek query names may carry .pdb or a chain suffix
    m = [q for q in Q if x == q or x.startswith(q + ".") or x.startswith(q + "_")]
    if len(m) != 1: raise SystemExit("ABORT: cannot map query name %s" % x)
    return m[0]
out = open(W + "/results/step29e_foldseek_summary.tsv", "w")
out.write("query\tgroup\tlength\tplddt\tptm\tdb\tn_hits\ttop_target\tbits\tprob\tevalue\tqtmscore\tttmscore\tqcov\tpass_qtm0.5_prob0.5\n")
for db in ("PDB", "CATH50"):
    best, n = {}, {}
    for l in open(W + "/results/step29e_foldseek_%s.m8" % db):
        f = l.rstrip("\n").split("\t")
        q = qid(f[0]); n[q] = n.get(q, 0) + 1
        if q not in best or float(f[7]) > float(best[q][7]): best[q] = f
    for q in Q:
        g = lab[q]; b = best.get(q)
        if b:
            ok = "yes" if float(b[10]) >= 0.5 and float(b[8]) >= 0.5 else "no"
            row = [b[1], b[7], b[8], b[6], b[10], b[11], b[13], ok]
        else:
            row = ["none"] + ["NA"] * 6 + ["no"]
        out.write("\t".join([q, g[1], g[2], g[4], g[5], db, str(n.get(q, 0))] + row) + "\n")
        print("%-32s %-10s %-6s n=%-4d %s" % (q, g[1], db, n.get(q, 0), " ".join(row)))
out.close()
print("wrote", W + "/results/step29e_foldseek_summary.tsv")
PYEOF
date
