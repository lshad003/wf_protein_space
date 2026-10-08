import gzip, collections, sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space/"
DROP = "UHM586.41010"
sup89, sup88 = set(), set()
def flush(rep, c):
    tot, ns = sum(c.values()), len(c)
    if tot >= 3 and ns >= 2: sup89.add(rep)
    c.pop(DROP, None)
    if sum(c.values()) >= 3 and len(c) >= 2: sup88.add(rep)
prev, c = None, collections.Counter()
for line in open(W + "results/clusters_95.tsv"):
    rep, mem = line.rstrip("\n").split("\t")
    if rep != prev:
        if prev is not None: flush(prev, c)
        prev, c = rep, collections.Counter()
    c[mem.split("__")[0]] += 1
flush(prev, c)
print("supported 89:", len(sup89), " supported 88:", len(sup88))
if len(sup89) != 2069453 or len(sup88) != 2067011: sys.exit("ABORT: support counts off")
with open(W + "results/supported_reps_95_n88.txt", "w") as o:
    for r in sorted(sup88): o.write(r + "\n")
fi = gzip.open(W + "results/count_matrix_primary.tsv.gz", "rt")
fo = gzip.open(W + "results/count_matrix_primary_n88.tsv.gz", "wt")
h = fi.readline().rstrip("\n").split("\t")
if h.count(DROP) != 1: sys.exit("ABORT: UHM586 column not found once")
k = [i for i, x in enumerate(h) if x != DROP]
fo.write("\t".join(h[i] for i in k) + "\n")
nin = nout = 0
for line in fi:
    f = line.rstrip("\n").split("\t"); nin += 1
    if f[0] in sup88:
        fo.write("\t".join(f[i] for i in k) + "\n"); nout += 1
fo.close()
print("matrix rows in:", nin, " out:", nout, " columns out:", len(k) - 1)
if nin != 2069453 or nout != 2067011: sys.exit("ABORT: matrix rows off (row IDs may not match rep IDs)")
m = [l for l in open(W + "metadata/samples_89_treatment.tsv") if not l.startswith(DROP + "\t")]
open(W + "metadata/samples_88_treatment.tsv", "w").writelines(m)
print("metadata rows (with header):", len(m))
print("DONE step25a")
