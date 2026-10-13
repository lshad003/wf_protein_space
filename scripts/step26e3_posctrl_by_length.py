# Step 26e3: positive-control own-family recovery by query length, and EU/GU query length distribution.
# Control: results/nmpfams2_posctrl_queries.faa (5,000 non-representative NMPFams2 metag members, step26e2; header
#   q<i>|<own family>) and results/nmpfams2_posctrl_hits.tsv (query,target,...; step26b3 settings: E <= 1e-3,
#   --min-seq-id 0.3, -c 0.7, --cov-mode 0, vs representatives). Own family = own family rep among the query's hits,
#   as in step26e2.
# Dark queries: results/unk_filtered_aplakidou.faa, class from results/fourway_classes_n88.tsv (col 2).
# Length = residues, terminal "*" removed (as step26e2). Bins: <100, 100-199, 200-399, >=400 aa.
# Read-only on every _n88 file.
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results"
OUT = R + "/nmpfams2_posctrl_by_length.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")

def lens(p):
    L = {}; k = None; s = []
    for l in open(p):
        if l[0] == ">":
            if k: L[k] = len("".join(s).rstrip("*"))
            k = l[1:].split()[0]; s = []
        else: s.append(l.strip())
    if k: L[k] = len("".join(s).rstrip("*"))
    return L

BINS = [("<100", 0, 100), ("100-199", 100, 200), ("200-399", 200, 400), (">=400", 400, 10 ** 9)]
def b(x): return next(n for n, lo, hi in BINS if lo <= x < hi)

ql = lens(R + "/nmpfams2_posctrl_queries.faa")
hits = collections.defaultdict(set)
for l in open(R + "/nmpfams2_posctrl_hits.tsv"):
    f = l.split("\t", 2); hits[f[0]].add(f[1])
own = {q for q in ql if q.split("|", 1)[1] in hits[q]}
print("control queries:", len(ql), "| own family recovered:", len(own))

cls = {}
for i, l in enumerate(open(R + "/fourway_classes_n88.tsv")):
    if i == 0: continue
    f = l.split("\t", 2); cls[f[0]] = f[1]
dl = lens(R + "/unk_filtered_aplakidou.faa")
print("dark queries:", len(dl))

def pct(v, p):
    v = sorted(v); return v[min(len(v) - 1, int(p * len(v)))]

rows = [("set", "bin", "n", "of", "pct", "own_family_recovered", "pct_recovered")]
for name, L, rec in [("posctrl", ql, own), ("EU", {k: v for k, v in dl.items() if cls[k] == "EU"}, None),
                     ("GU", {k: v for k, v in dl.items() if cls[k] == "GU"}, None)]:
    c = collections.Counter(b(v) for v in L.values())
    for n, _, _ in BINS + [("all", 0, 0)]:
        ks = list(L) if n == "all" else [k for k in L if b(L[k]) == n]
        r = sum(k in rec for k in ks) if rec is not None else None
        rows.append((name, n, len(ks), len(L), "%.2f" % (100.0 * len(ks) / len(L)),
                     "NA" if r is None else r, "NA" if r is None else "%.2f" % (100.0 * r / len(ks)) if ks else "NA"))
    v = list(L.values())
    for q, p in (("min", 0.0), ("p10", 0.1), ("p25", 0.25), ("median", 0.5), ("p75", 0.75), ("p90", 0.9)):
        rows.append((name, q, pct(v, p), "NA", "NA", "NA", "NA"))
    rows.append((name, "max", max(v), "NA", "NA", "NA", "NA"))
    rows.append((name, "mean", "%.1f" % (sum(v) / float(len(v))), "NA", "NA", "NA", "NA"))
with open(OUT, "w") as fh:
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
for r in rows: print("\t".join(map(str, r)))
print("wrote", OUT)
