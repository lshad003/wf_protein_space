#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 128gb
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step2c_diag.%j.log
module load mmseqs2/13-45111
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
V1=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
TMP=$SCRATCH/step2c_$SLURM_JOB_ID
mkdir -p $TMP
CPU=${SLURM_CPUS_ON_NODE:-32}

echo "node: $(hostname)"
echo "cpu flags of interest: $(lscpu | grep -o 'avx512[a-z]*\|avx2\|sse4_1' | sort -u | tr '\n' ' ')"
echo

QRY=$WPS/catalog/subsample_200k.faa
echo "query seqs: $(grep -c '^>' $QRY)"

mmseqs easy-search $QRY $V1/db/LsFMGC_AA_95_rep.fasta \
  $WPS/results/step2c_hits.m8 $TMP \
  --min-seq-id 0.95 -c 0.8 --cov-mode 1 \
  --max-accept 1 -e 1e-5 -s 4 --threads $CPU

if [ ! -s $WPS/results/step2c_hits.m8 ]; then
  echo "SEARCH FAILED, no output file"
  exit 1
fi

echo "== RESULT =="
Q=$(grep -c '^>' $QRY)
H=$(cut -f1 $WPS/results/step2c_hits.m8 | sort -u | wc -l)
echo "queries: $Q ; with hit: $H"
$PY -c "print('OVERALL mapping rate: %.2f%%' % (100.0*$H/$Q))"

echo "== per-year =="
grep '^>' $QRY | sed 's/>//; s/__.*//' > $TMP/q_stems.txt
cut -f1 $WPS/results/step2c_hits.m8 | sort -u | sed 's/__.*//' > $TMP/h_stems.txt
$PY << PYEOF
from collections import Counter
q = Counter(open("$TMP/q_stems.txt").read().split())
h = Counter(open("$TMP/h_stems.txt").read().split())
wf23 = set(open("$WPS/metadata/wf23_stems.txt").read().split())
tot = {"WF23":[0,0], "WF24":[0,0]}
rows = []
for s in sorted(q):
    yr = "WF23" if s in wf23 else "WF24"
    tot[yr][0] += q[s]; tot[yr][1] += h.get(s,0)
    rows.append((s, yr, q[s], h.get(s,0), 100.0*h.get(s,0)/q[s]))
for yr in ("WF23","WF24"):
    a,b = tot[yr]
    print("%s: %d queries, %d hits, %.2f%%" % (yr, a, b, 100.0*b/a))
print()
print("lowest 5 samples by mapping rate:")
for r in sorted(rows, key=lambda x: x[4])[:5]:
    print("  %s %s %d/%d = %.1f%%" % (r[0], r[1], r[3], r[2], r[4]))
PYEOF
rm -rf $TMP
