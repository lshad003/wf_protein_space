#!/bin/bash
#SBATCH -J step34n_antifam
#SBATCH -p short -c 8 --mem 16G --time=02:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34n_antifam_eu.%j.log
# Step 34n: AntiFam screen of the 41,114 strict_final EU representatives (results/step34a_rep_origin_v2_final.tsv,
#   class4 EU, strict_final 1). AntiFam release 8.0 (db/antifam, 278 families, 2023-10-12; the copy step14b used),
#   HMMER hmmsearch --cut_ga. Sequences from catalog/db/LsPS_AA_95_rep.fasta, terminal "*" stripped (count printed).
# Outputs: results/step34n/eu_reps.faa, results/step34n_antifam_eu.tblout, results/step34n_antifam_eu.domtblout.
# Exits with the main commands' status (CLAUDE.md rule); refuses to overwrite.
W=/bigdata/stajichlab/lshad003/wf_protein_space
AF=$W/db/antifam; O=$W/results; D=$O/step34n
for f in $O/step34n_antifam_eu.tblout $O/step34n_antifam_eu.domtblout $D/eu_reps.faa; do
  if [ -e $f ]; then echo "REFUSING: $f exists"; exit 1; fi
done
[ -s $AF/AntiFam.hmm ] || { echo "ABORT: $AF/AntiFam.hmm missing"; exit 1; }
mkdir -p $D
hostname; date
echo "AntiFam: $(head -1 $AF/version); families $(grep -c '^NAME' $AF/AntiFam.hmm)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - <<'PYEOF'
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
eu = set()
with open(W + "/results/step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, isf = h.index("class4"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[isf] == "1" and f[ic] == "EU": eu.add(f[0])
if len(eu) != 41114: raise SystemExit("ABORT: %d EU reps" % len(eu))
seq = {}; cur = None
for l in open(W + "/catalog/db/LsPS_AA_95_rep.fasta"):
    if l[0] == ">":
        cur = l[1:].split()[0]; cur = cur if cur in eu else None
        if cur: seq[cur] = []
    elif cur: seq[cur].append(l.strip())
star = internal = 0
with open(W + "/results/step34n/eu_reps.faa", "w") as o:
    for r in sorted(seq):
        s = "".join(seq[r])
        if s.endswith("*"): s = s[:-1]; star += 1
        if "*" in s: internal += 1
        o.write(">%s\n%s\n" % (r, s))
print("EU reps listed %d, sequences written %d, terminal * stripped %d, internal * %d" % (len(eu), len(seq), star, internal))
if len(seq) != len(eu): raise SystemExit("ABORT: %d EU reps missing from LsPS_AA_95_rep.fasta" % (len(eu) - len(seq)))
PYEOF
rc=$?
if [ $rc -ne 0 ]; then echo "ABORT: extraction exit $rc"; exit 1; fi
module load hmmer/3.4
hmmsearch -h | head -2 | tail -1
hmmsearch --cut_ga --cpu $SLURM_CPUS_PER_TASK --tblout $O/step34n_antifam_eu.tblout --domtblout $O/step34n_antifam_eu.domtblout \
  $AF/AntiFam.hmm $D/eu_reps.faa > /dev/null
rc=$?
echo "hmmsearch exit $rc"
if [ $rc -eq 0 ]; then
  echo "EU reps with an AntiFam hit (--cut_ga): $(grep -v '^#' $O/step34n_antifam_eu.tblout | awk '{print $1}' | sort -u | wc -l) of 41114"
  echo "top families:"; grep -v '^#' $O/step34n_antifam_eu.tblout | awk '{print $3}' | sort | uniq -c | sort -rn | head -10
fi
date
exit $rc
