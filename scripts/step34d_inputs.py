# Step 34d inputs: ProstT5/Foldseek query sets for the strict_final recurrent EU families (method as step28h).
# Families: results/step34c2c_strict_v2_families.tsv, recurrent = 1 and excluded = 0 (4,509). Query = the family
#   representative (Step 34c 30% cluster rep) if it is in strict_final, else the longest strict_final EU rep of the
#   family (ties by ID); strict_final from results/step34a_rep_origin_v2_final.tsv; family reps from
#   results/step34c_strict_eu_clusters30.tsv. Sequences from catalog/db/LsPS_AA_95_rep.fasta, terminal "*" stripped
#   (aborts on an internal "*").
# Controls: (a) 4,509 strict_final K reps drawn at random (seed 34) within 10-aa length bins to match the query
#   length distribution; where a bin has too few K reps, all are taken and the rest borrowed from the nearest bins
#   (upper first; counts printed); (b) every query reversed, headers "rev_".
# Outputs (refuse to overwrite): results/step34d_queries.faa, results/step34d_known_lenmatched.faa,
#   results/step34d_decoys_rev.faa, results/step34d_query_map.tsv
import os, random, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OQ, OK, OD, OM = R + "step34d_queries.faa", R + "step34d_known_lenmatched.faa", R + "step34d_decoys_rev.faa", R + "step34d_query_map.tsv"
for p in (OQ, OK, OD, OM):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
EU, K = set(), []
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, ifn = h.index("class4"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ifn] == "1":
            if f[ic] == "EU": EU.add(f[0])
            elif f[ic] == "K": K.append(f[0])
if len(EU) != 41114: raise SystemExit("ABORT: strict_final EU %d" % len(EU))
fams = []
with open(R + "step34c2c_strict_v2_families.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ir, ix = h.index("recurrent"), h.index("excluded")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ir] == "1" and f[ix] == "0": fams.append(f[0])
if len(fams) != 4509: raise SystemExit("ABORT: %d recurrent families" % len(fams))
reps_of = collections.defaultdict(list)
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t")
    if b in EU: reps_of[a].append(b)
need = set(EU) | set(K)
seq = {}; hd = None; buf = []
with open(W + "/catalog/db/LsPS_AA_95_rep.fasta") as fh:
    for l in fh:
        if l[0] == ">":
            if hd in need: seq[hd] = "".join(buf)
            hd = l[1:].split()[0]; buf = []
        elif hd in need: buf.append(l.strip())
    if hd in need: seq[hd] = "".join(buf)
if len(seq) != len(need): raise SystemExit("ABORT: %d sequences missing" % (len(need) - len(seq)))
star = sum(1 for v in seq.values() if v.endswith("*"))
seq = {k: v.rstrip("*") for k, v in seq.items()}
bad = [k for k, v in seq.items() if "*" in v]
if bad: raise SystemExit("ABORT: internal stop in %d sequences" % len(bad))
print("sequences read %d (EU %d + K %d); terminal '*' stripped from %d" % (len(seq), len(EU), len(K), star))
q = []; fb = 0
for f in fams:
    if f in EU: q.append((f, f, "family_rep"))
    else:
        r = sorted(reps_of[f], key=lambda x: (-len(seq[x]), x))[0]; q.append((f, r, "longest_strict_final")); fb += 1
print("queries %d | family representative %d | fallback to longest strict_final EU rep %d" % (len(q), len(q) - fb, fb))
L = [len(seq[r]) for _, r, _ in q]
bn = lambda n: n // 10
needb = collections.Counter(bn(n) for n in L)
pool = collections.defaultdict(list)
for k in sorted(K): pool[bn(len(seq[k]))].append(k)
rng = random.Random(34); kp = []; borrowed = 0; used = set()
for b in sorted(needb):
    take = rng.sample(pool[b], min(needb[b], len(pool[b]))); kp += take; used |= set(take)
    short = needb[b] - len(take); d = 1
    while short > 0:   # bin too small: borrow from the nearest bins, upper first
        for nb in (b + d, b - d):
            av = [k for k in pool.get(nb, []) if k not in used]
            if short > 0 and av:
                t = rng.sample(av, min(short, len(av))); kp += t; used |= set(t); short -= len(t); borrowed += len(t)
                print("  10-aa bin %d-%d aa short: borrowed %d from bin %d-%d aa" % (b * 10, b * 10 + 9, len(t), nb * 10, nb * 10 + 9))
        d += 1
        if d > 100: raise SystemExit("ABORT: cannot fill bin %d" % b)
print("K drawn %d; borrowed from adjacent 10-aa bins %d" % (len(kp), borrowed))
with open(OQ, "w") as o:
    for _, r, _ in q: o.write(">%s\n%s\n" % (r, seq[r]))
with open(OK, "w") as o:
    for k in kp: o.write(">%s\n%s\n" % (k, seq[k]))
with open(OD, "w") as o:
    for _, r, _ in q: o.write(">rev_%s\n%s\n" % (r, seq[r][::-1]))
with open(OM, "w") as o:
    o.write("family_rep\tquery\tsource\tlength_aa\n")
    for f, r, s in q: o.write("%s\t%s\t%s\t%d\n" % (f, r, s, len(seq[r])))
def d(lab, X):
    X = sorted(X); print("%s n=%d min %d median %d mean %.1f max %d | <100 %d, 100-199 %d, 200-399 %d, >=400 %d" % (lab, len(X), X[0], X[len(X) // 2],
          sum(X) / len(X), X[-1], sum(x < 100 for x in X), sum(100 <= x < 200 for x in X), sum(200 <= x < 400 for x in X), sum(x >= 400 for x in X)))
d("queries", L); d("K lenmatched", [len(seq[k]) for k in kp]); d("decoys", L)
print("wrote", OQ, OK, OD, OM)
