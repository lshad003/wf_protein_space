# Step 34d2 inputs: fix the Step 34d K control. results/step34d_known_lenmatched.faa has 4,509 records but 225 IDs
# twice (40-49 aa; the 40-49 aa bin was drawn without excluding the 246 reps borrowed for 30-39 aa). The second copy of
# each duplicate is replaced by a fresh strict_final K rep (seed 3442), drawn excluding every rep already in the K set:
# all available 40-49 aa reps first (27), the rest from the nearest upper 10-aa bins (50-59 aa, then up), the
# nearest-bin rule of step34d_inputs.py (option 1, agreed 2026-10-16). Sequences from catalog/db/LsPS_AA_95_rep.fasta, terminal "*" stripped.
# Outputs (refuse to overwrite; old files kept): results/step34d2_k_replacements.faa (225),
#   results/step34d_known_lenmatched_fixed.faa (4,509 distinct, same order, replacements in place of second copies).
import os, random, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
O1, O2 = R + "step34d2_k_replacements.faa", R + "step34d_known_lenmatched_fixed.faa"
for p in (O1, O2):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
rec = []; h = None
for l in open(R + "step34d_known_lenmatched.faa"):
    if l[0] == ">": h = l[1:].strip()
    else: rec.append((h, l.strip()))
c = collections.Counter(k for k, _ in rec); dup = {k for k, v in c.items() if v > 1}
print("K records %d | distinct %d | duplicated IDs %d (lengths %d-%d)" % (len(rec), len(c), len(dup),
      min(len(s) for k, s in rec if k in dup), max(len(s) for k, s in rec if k in dup)))
if len(rec) != 4509 or len(dup) != 225 or max(c.values()) != 2: raise SystemExit("ABORT: unexpected duplicate structure")
K = set()
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    hd = next(fh).rstrip("\n").split("\t"); ic, ifn = hd.index("class4"), hd.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ifn] == "1" and f[ic] == "K": K.add(f[0])
seq = {}; hd = None; buf = []
with open(W + "/catalog/db/LsPS_AA_95_rep.fasta") as fh:
    for l in fh:
        if l[0] == ">":
            if hd in K: seq[hd] = "".join(buf).rstrip("*")
            hd = l[1:].split()[0]; buf = []
        elif hd in K: buf.append(l.strip())
    if hd in K: seq[hd] = "".join(buf).rstrip("*")
used = set(c)
rng = random.Random(3442); new = []; lo = 40
while len(new) < 225:
    pool = sorted(k for k, s in seq.items() if lo <= len(s) <= lo + 9 and k not in used)
    t = rng.sample(pool, min(225 - len(new), len(pool))); new += t; used |= set(t)
    print("strict_final K reps %d-%d aa not in the K set: %d | drawn %d" % (lo, lo + 9, len(pool), len(t)))
    lo += 10
    if lo >= 100: raise SystemExit("ABORT: replacements would leave the <100 aa bin")
if any("*" in seq[k] for k in new): raise SystemExit("ABORT: internal stop")
it = iter(new); seen = set(); out = []
for k, s in rec:
    if k in dup and k in seen: n = next(it); out.append((n, seq[n]))
    else: out.append((k, s)); seen.add(k)
if len({k for k, _ in out}) != 4509: raise SystemExit("ABORT: fixed set not 4,509 distinct")
with open(O1, "w") as o:
    for k in new: o.write(">%s\n%s\n" % (k, seq[k]))
with open(O2, "w") as o:
    for k, s in out: o.write(">%s\n%s\n" % (k, s))
q = [len(l.strip()) for l in open(R + "step34d_queries.faa") if l[0] != ">"]
b = lambda X: [sum(x < 100 for x in X), sum(100 <= x < 200 for x in X), sum(200 <= x < 400 for x in X), sum(x >= 400 for x in X)]
X = [len(s) for _, s in out]
print("fixed K: n %d distinct %d median %d | bins %s | queries bins %s" % (len(out), len({k for k, _ in out}), sorted(X)[len(X) // 2], b(X), b(q)))
print("10-aa bin counts equal to queries in every bin:", collections.Counter(x // 10 for x in X) == collections.Counter(x // 10 for x in q))
print("wrote", O1, "and", O2)
