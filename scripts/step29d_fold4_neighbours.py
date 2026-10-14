import glob, collections, re
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
NAME = {"UHM125.35769__k141_151130_2": "DivIVA-like", "UHM636.41030__k141_332084_5": "Cse2/CasB-like",
        "UHM145.35772__k141_106404_54": "HTH-like (weak)", "UHM145.35771__k141_106979_4": "helical repeat-like (weak)"}
OFF = (-3, -2, -1, 1, 2, 3)
mem = {q: [] for q in NAME}
for l in open(W + "/results/unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in mem: mem[r].append(m)
def nb(g, off):
    stem, rest = g.split("__", 1); i = rest.rfind("_")
    try: n = int(rest[i + 1:])
    except ValueError: return None
    return "%s__%s_%d" % (stem, rest[:i], n + off)
need = {nb(m, o) for ms in mem.values() for m in ms for o in OFF} - {None}
rep = {}
for l in open(W + "/results/clusters_95.tsv"):
    r, m = l.rstrip("\n").split("\t")[:2]
    if m in need: rep[m] = r
reps = set(rep.values()); dom = collections.defaultdict(set)
for f in glob.glob(W + "/results/pfam_hs/*.tblout"):
    for l in open(f):
        if l[0] == "#": continue
        p = l.split()
        if p[0] in reps: dom[p[0]].add(p[2])
print("members %s | neighbour slots %d, genes found %d, with Pfam %d" % ({NAME[q]: len(v) for q, v in mem.items()}, len(need), len(rep), sum(1 for n in rep if dom.get(rep[n]))))
out = open(W + "/results/step29d_fold4_neighbours.tsv", "w"); out.write("family_rep\tlabel\tmember\toffset\tneighbour\trep95\tpfam\n")
for q in NAME:
    w3, w1, anyg, cas = collections.Counter(), collections.Counter(), 0, 0
    for m in mem[q]:
        s3, s1, has = set(), set(), False
        for o in OFF:
            n = nb(m, o)
            if n not in rep: continue
            has = True; d = dom.get(rep[n], set())
            s3 |= d; s1 |= d if abs(o) == 1 else set()
            out.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\n" % (q, NAME[q], m, o, n, rep[n], ",".join(sorted(d)) or "-"))
        anyg += has; w3.update(s3); w1.update(s1)
        cas += any(re.search(r"cas|crispr|cse", x, re.I) for x in s3)
    print("\n== %s  %s | %d members, %d with a neighbour gene | members with a CRISPR/Cas-named domain within 3 genes: %d" % (q, NAME[q], len(mem[q]), anyg, cas))
    print("  within 1 gene :", ", ".join("%s %d" % x for x in w1.most_common(10)) or "none")
    print("  within 3 genes:", ", ".join("%s %d" % x for x in w3.most_common(15)) or "none")
out.close()
