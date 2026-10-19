# Step 34j A2d: whole-window context similarity, a second conservation score next to the top-neighbour fraction (which
#   saturates: control median 1.000 under A2, A2b and A2c, logs/step34j_a2*.log).
# Over the containment-collapsed backgrounds of A2c (results/step34j/context_nr50c_seeds.tsv): each seed's set of
#   neighbour 50% clusters within +-5 (own 50% cluster excluded, as A2); score = mean pairwise Jaccard over all seed pairs
#   (units with >= 2 backgrounds). Percentile: mid-rank among controls in the same background-count bin (2, 3-4, 5-9, >=10).
# Output (refuses to overwrite): results/step34j/context_jaccard.tsv
import os, collections, itertools
D = "/bigdata/stajichlab/lshad003/wf_protein_space/results/step34j/"; O = D + "context_jaccard.tsv"
if os.path.exists(O): raise SystemExit("REFUSING: " + O + " exists")
seeds = collections.defaultdict(list)
with open(D + "context_nr50c_seeds.tsv") as fh:
    next(fh)
    for l in fh:
        st, f, m, n = l.rstrip("\n").split("\t"); seeds[(st, f)].append(m)
want = {m for v in seeds.values() for m in v}
own = {}
with open("/bigdata/stajichlab/lshad003/wf_protein_space/results/clusters_50_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m in want: own[m] = r
nb = collections.defaultdict(set)
with open(D + "context_neighbours.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[2] in want and f[8] != own.get(f[2]): nb[(f[0], f[1], f[2])].add(f[8])
def bn(n): return "1" if n < 2 else "2" if n == 2 else "3-4" if n <= 4 else "5-9" if n <= 9 else ">=10"
res = {}
for k, ms in seeds.items():
    if len(ms) < 2: res[k] = (len(ms), None); continue
    sets = [nb.get((k[0], k[1], m), set()) for m in ms]
    js = [len(a & b) / len(a | b) if (a | b) else 0.0 for a, b in itertools.combinations(sets, 2)]
    res[k] = (len(ms), sum(js) / len(js))
ctl = collections.defaultdict(list)
for (st, f), (n, j) in res.items():
    if st == "control" and j is not None: ctl[bn(n)].append(j)
def pct(j, b):
    v = ctl.get(b)
    return "%.1f" % (100.0 * (sum(x < j for x in v) + 0.5 * sum(x == j for x in v)) / len(v)) if v else "NA"
with open(O, "w") as o:
    o.write("set\tfamily\tn_backgrounds\tbackground_bin\tmean_pairwise_jaccard\tjaccard_percentile_vs_controls\n")
    for (st, f), (n, j) in sorted(res.items()):
        o.write("%s\t%s\t%d\t%s\t%s\t%s\n" % (st, f, n, bn(n), "%.3f" % j if j is not None else "NA", pct(j, bn(n)) if j is not None else "NA"))
for st in ("EU", "control"):
    for b in ("2", "3-4", "5-9", ">=10"):
        v = sorted(j for (s, f), (n, j) in res.items() if s == st and j is not None and bn(n) == b)
        if v: print("%-7s bin %-4s n %4d | quartiles %.3f %.3f %.3f | 90th %.3f" % (st, b, len(v), v[len(v) // 4], v[len(v) // 2], v[3 * len(v) // 4], v[int(0.9 * (len(v) - 1))]))
print("wrote", O)
