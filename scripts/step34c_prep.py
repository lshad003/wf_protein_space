# Step 34c prep: FASTA of the strict EU reps for 30% clustering.
# Strict EU = class4 EU and origin Tiara-prokaryotic in results/step34a_rep_origin.tsv (Step 34a, 41,392).
# Sequences from catalog/db/LsPS_AA_95_rep.fasta with the terminal "*" stripped (Step 28h lesson).
# Usage: python3 step34c_prep.py <out.faa>; aborts unless exactly 41,392 sequences are written.
import sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
WANT = 41392
want = set()
with open(W + "/results/step34a_rep_origin.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[1] == "EU" and f[2] == "Tiara-prokaryotic": want.add(f[0])
print("strict EU reps in step34a_rep_origin.tsv:", len(want), flush=True)
if len(want) != WANT: raise SystemExit("ABORT: expected %d strict EU reps" % WANT)
n = 0; star = 0; inner = 0; dup = 0; got = set()
out = open(sys.argv[1], "w")
def emit(h, s):
    global n, star, inner, dup
    if h in got: dup += 1; return
    if s.endswith("*"): star += 1; s = s.rstrip("*")
    if "*" in s: inner += 1
    got.add(h); n += 1; out.write(">%s\n%s\n" % (h, s))
h = None; s = []
with open(W + "/catalog/db/LsPS_AA_95_rep.fasta") as fh:
    for l in fh:
        if l[0] == ">":
            if h is not None and h in want: emit(h, "".join(s))
            h = l[1:].split()[0]; s = []
        elif h in want: s.append(l.strip())
if h is not None and h in want: emit(h, "".join(s))
out.close()
print("written:", n, "| terminal * stripped:", star, "| internal *:", inner, "| duplicate IDs:", dup, flush=True)
if n != WANT or dup or inner: raise SystemExit("ABORT: written %d, duplicates %d, internal * %d" % (n, dup, inner))
