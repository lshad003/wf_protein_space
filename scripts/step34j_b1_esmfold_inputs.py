# Step 34j B1 inputs: ESMFold query set.
# EU: the Step 34d query (one rep per family, results/step34d_queries.faa) of each of the 4,250 member-validated
#   recurrent families (step34c2c_strict_v2_families.tsv recurrent == 1, excluded == 0; every 95% cluster of the family
#   has class_any == EU in results/step34h_member_classes.tsv).
# Controls: all 4,509 known length-matched K reps in results/step34d_known_lenmatched_fixed.faa (none has an
#   existing model in results/esmfold_pdb; checked here).
# Output (refuses to overwrite): results/step34j/esmfold_input.faa, results/step34j/esmfold_input_labels.tsv,
#   results/step34j/esmfold_chunks/C.<n>.faa (N chunks, length-sorted round robin so tasks carry similar load).
import os, sys, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; O = R + "step34j/"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for p in (O + "esmfold_input.faa", O + "esmfold_input_labels.tsv", O + "esmfold_chunks"):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def fa(p):
    d, h = {}, None
    for l in open(p):
        if l[0] == ">": h = l[1:].split()[0]; d[h] = []
        else: d[h].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}
fams = set()
with open(R + "step34c2c_strict_v2_families.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ir, ix = h.index("recurrent"), h.index("excluded")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ir] == "1" and f[ix] == "0": fams.add(f[0])
cls = collections.defaultdict(list)
with open(R + "step34h_member_classes.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ifr, ic = h.index("family_rep"), h.index("class_any")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ifr] in fams: cls[f[ifr]].append(f[ic])
mv = {f for f in fams if cls[f] and all(c == "EU" for c in cls[f])}
print("recurrent families %d, member-validated %d" % (len(fams), len(mv)))
if len(fams) != 4509 or len(mv) != 4250: raise SystemExit("ABORT: expected 4,509 and 4,250")
qmap = {}
with open(R + "step34d_query_map.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); qmap[f[0]] = f[1]
qs, ks = fa(R + "step34d_queries.faa"), fa(R + "step34d_known_lenmatched_fixed.faa")
if len(ks) != 4509: raise SystemExit("ABORT: controls %d" % len(ks))
have = {x[:-4] for x in os.listdir(R + "esmfold_pdb") if x.endswith(".pdb")}
print("controls with an existing model in results/esmfold_pdb: %d" % len(have & set(ks)))
rows = [(qmap[f], qs[qmap[f]], "EU", f) for f in sorted(mv)] + [(k, s, "control", "-") for k, s in sorted(ks.items()) if k not in have]
ids = [r[0] for r in rows]
if len(ids) != len(set(ids)): raise SystemExit("ABORT: duplicate IDs between sets")
bad = [r[0] for r in rows if "*" in r[1] or not r[1]]
if bad: raise SystemExit("ABORT: stop or empty in %d sequences" % len(bad))
os.makedirs(O + "esmfold_chunks")
with open(O + "esmfold_input.faa", "w") as o, open(O + "esmfold_input_labels.tsv", "w") as t:
    t.write("id\tset\tfamily_rep\tlength_aa\tchunk\n")
    rows.sort(key=lambda r: (-len(r[1]), r[0]))
    ch = {}
    for i, r in enumerate(rows):
        c = i % N; ch.setdefault(c, []).append(r)
        o.write(">%s\n%s\n" % (r[0], r[1])); t.write("%s\t%s\t%s\t%d\t%d\n" % (r[0], r[2], r[3], len(r[1]), c))
for c, rs in ch.items():
    with open(O + "esmfold_chunks/C.%d.faa" % c, "w") as o:
        for r in sorted(rs, key=lambda r: len(r[1])): o.write(">%s\n%s\n" % (r[0], r[1]))
L = sorted(len(r[1]) for r in rows)
print("queries %d (EU %d, control %d); length min %d median %d max %d; chunks %d" % (len(rows), sum(r[2] == "EU" for r in rows),
      sum(r[2] == "control" for r in rows), L[0], L[len(L) // 2], L[-1], N))
