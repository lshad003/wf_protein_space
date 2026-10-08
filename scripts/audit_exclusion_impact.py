import collections, sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space/"
X = ["UHM586.41010", "UHM590.41012"]
SC = {"UHM586": {X[0]}, "UHM590": {X[1]}, "both": set(X)}
prot = collections.Counter(); contain = collections.Counter(); only = collections.Counter()
lost = collections.Counter(); emptied = collections.Counter()
repfrom = collections.Counter(); repfrom_sup = collections.Counter()
n = {"cl": 0, "sup": 0}
def flush(rep, c):
    n["cl"] += 1
    tot, ns = sum(c.values()), len(c)
    s = tot >= 3 and ns >= 2
    n["sup"] += s
    rs = rep.split("__")[0]
    if rs in X:
        repfrom[rs] += 1; repfrom_sup[rs] += s
    for st in X:
        if st in c:
            contain[st] += 1; only[st] += (ns == 1)
    for k, rem in SC.items():
        t2 = tot - sum(c.get(r, 0) for r in rem)
        n2 = ns - sum(1 for r in rem if r in c)
        if t2 == 0: emptied[k] += 1
        if s and not (t2 >= 3 and n2 >= 2): lost[k] += 1
prev, c = None, collections.Counter()
with open(W + "results/clusters_95.tsv") as f:
    for line in f:
        rep, mem = line.rstrip("\n").split("\t")
        if rep != prev:
            if prev is not None: flush(prev, c)
            prev, c = rep, collections.Counter()
        st = mem.split("__")[0]; c[st] += 1; prot[st] += 1
flush(prev, c)
print("stems:", len(prot), "clusters:", n["cl"], "supported:", n["sup"])
if (len(prot), n["cl"], n["sup"]) != (89, 6182117, 2069453):
    sys.exit("ABORT: totals do not reproduce known catalog; clusters may not be contiguous by rep")
meta = {l.split("\t")[0]: l.split("\t")[1] for l in open(W + "metadata/samples_89_treatment.tsv")}
wf24 = sorted(prot[s] for s in prot if meta.get(s) == "WF24")
print("WF24 proteins per metagenome: min", wf24[0], "median", wf24[len(wf24)//2], "max", wf24[-1])
for st in X:
    rank = sum(1 for v in wf24 if v < prot[st]) + 1
    print(st, "proteins", prot[st], "WF24 rank", rank, "of", len(wf24), "| clusters containing", contain[st],
          "| clusters with only this stem", only[st], "| reps from it", repfrom[st], "supported reps", repfrom_sup[st])
for k in SC:
    print("remove", k, "| clusters emptied", emptied[k], "| supported clusters losing support", lost[k],
          "(%.4f%% of 2,069,453)" % (100 * lost[k] / 2069453))
print("DONE audit_exclusion_impact")
