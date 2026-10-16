# Step 33g summary: D. melanogaster (dm6) status of the step33c contig sets from results/step33g_paf (minimap2 2.30
# -x asm10 --secondary=no -c). Per contig over all PAF rows (every row must carry cg:Z): coverage = merged query
# interval length / contig length; identity = sum(col 10, base-level matches) / sum(col 11, block length).
# Fly = coverage >= 0.50 and identity >= 0.95 (main); identity >= 0.90 sensitivity (rule as step33c2).
# Reports: A (Tiara-eukaryotic >= 1 kb), B (contigs with an EU supported rep; as defined incl. the 405 EU reps on
# UHM586.41010, and excluding them; per contig and per EU rep), C (negative control, information only), F (contigs with
# a member of the 15,580 families) and families with >= 50% of members on fly contigs; Cse2 family (expected 0);
# overlap with the step33c2 frog calls (results/step33c_contig_host_v2.tsv).
# Outputs (refuse to overwrite): results/step33g_contig_fly.tsv, results/step33g_family_fly.tsv
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; L = R + "step33c_lists/"; P = R + "step33g_paf/"
OUT = [R + "step33g_contig_fly.tsv", R + "step33g_family_fly.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
X = "UHM586.41010"; CSE2 = "UHM636.41030__k141_332084_5"
stems = [l.strip() for l in open(L + "stems.txt") if l.strip()]
if len(stems) != 89 or stems[-1] != X: raise SystemExit("ABORT: stems.txt is not the 88 + " + X)
sets = {}
for s in stems:
    for l in open(L + s + ".tsv"): c, v = l.rstrip("\n").split("\t"); sets[(s, c)] = v
iv, mt, bl, ql = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
for s in stems:
    p = P + s + ".fly.paf"
    if not os.path.exists(p): raise SystemExit("ABORT: missing " + p)
    for l in open(p):
        f = l.split("\t"); k = (s, f[0])
        if k not in sets: raise SystemExit("ABORT: unlisted contig in " + p)
        if "\tcg:Z:" not in l: raise SystemExit("ABORT: PAF row without cg:Z in " + p)
        ql[k] = int(f[1]); iv[k].append((int(f[2]), int(f[3]))); mt[k] += int(f[9]); bl[k] += int(f[10])
st = {}
for k, v in iv.items():
    v.sort(); c, cs, ce = 0, -1, -1
    for a, b in v:
        if a > ce: c += max(0, ce - cs); cs, ce = a, b
        else: ce = max(ce, b)
    st[k] = ((c + ce - cs) / ql[k], mt[k] / bl[k], ql[k])
def call(k, t): v = st.get(k); return bool(v and v[0] >= 0.5 and v[1] >= t)
frog = {}
for i, l in enumerate(open(R + "step33c_contig_host_v2.tsv")):
    if i: f = l.rstrip("\n").split("\t"); frog[(f[0], f[1])] = f[8] == "1"
by = {x: [k for k, v in sets.items() if x in v] for x in "ABCF"}
by["B_excl"] = [k for k in by["B"] if k[0] != X]
print("fly = coverage >= 50% and identity >= 95% (CIGAR); sensitivity identity >= 90%")
for x, lab in (("A", "Tiara-eukaryotic >= 1 kb"), ("B", "contigs with an EU supported rep (as defined)"), ("B_excl", "B excluding contigs on " + X),
               ("C", "negative control: 2% Tiara-prokaryotic >= 3 kb (info)"), ("F", "contigs with a member of the 15,580 families")):
    ks = by[x]; n = len(ks); a = sum(call(k, .95) for k in ks); b = sum(call(k, .90) for k in ks)
    print("  %-6s %-55s n=%d | fly95 %d (%.2f%%) | fly90 %d (%.2f%%) | any alignment %d | fly95 and frog95 %d"
          % (x, lab, n, a, 100.0 * a / n, b, 100.0 * b / n, sum(1 for k in ks if k in st), sum(1 for k in ks if call(k, .95) and frog.get(k))))
with open(OUT[0], "w") as o:
    o.write("stem\tcontig\tsets\tlength\tfly_cov\tfly_id\tfly95\tfly90\n")
    for k, v in sorted(sets.items()):
        h = st.get(k); o.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n" % (k[0], k[1], v, h[2] if h else "NA", "%.4f" % h[0] if h else "0",
                               "%.4f" % h[1] if h else "NA", call(k, .95), call(k, .90)))
def contig(g): s, rest = g.split("__", 1); return s, rest[:rest.rfind("_")]
e = [0, 0, 0]; ex = [0, 0, 0]
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        r, c = l.split("\t", 2)[:2]
        if c != "EU": continue
        k = contig(r)
        for t in ((e,) if k[0] == X else (e, ex)): t[0] += 1; t[1] += call(k, .95); t[2] += call(k, .90)
print("  per EU supported rep: as defined %d/%d (%.3f%%) fly95, %d fly90 | excluding %s %d/%d (%.3f%%) fly95, %d fly90"
      % (e[1], e[0], 100.0 * e[1] / e[0], e[2], X, ex[1], ex[0], 100.0 * ex[1] / ex[0], ex[2]))
fam = set()
with open(R + "unk_direct30_eu_recurrence.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ia, ic = h.index("n_animals"), h.index("n_cohorts")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if int(f[ia]) >= 10 and int(f[ic]) >= 2: fam.add(f[0])
if len(fam) != 15580: raise SystemExit("ABORT: %d families" % len(fam))
fm = {f: [0, 0, 0] for f in fam}
for l in open(R + "unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in fm: k = contig(m); v = fm[r]; v[0] += 1; v[1] += call(k, .95); v[2] += call(k, .90)
m95 = sum(1 for v in fm.values() if 2 * v[1] >= v[0]); m90 = sum(1 for v in fm.values() if 2 * v[2] >= v[0]); any95 = sum(1 for v in fm.values() if v[1])
with open(OUT[1], "w") as o:
    o.write("family_rep\tn_members\tn_fly95\tn_fly90\tmajority_fly95\tmajority_fly90\n")
    for f in sorted(fm): v = fm[f]; o.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (f, v[0], v[1], v[2], 2 * v[1] >= v[0], 2 * v[2] >= v[0]))
print("\n15,580 families (majority fly = fly members >= 50%% of members): fly95 %d (%.2f%%), fly90 %d (%.2f%%); with >= 1 fly95 member %d"
      % (m95, 100.0 * m95 / 15580, m90, 100.0 * m90 / 15580, any95))
v = fm[CSE2]; print("Cse2 family %s: %d members, fly95 %d, fly90 %d (expected 0): %s" % (CSE2, v[0], v[1], v[2], "PASS" if v[1] == v[2] == 0 else "FAIL"))
for p in OUT: print("wrote", p)
