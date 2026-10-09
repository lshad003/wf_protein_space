# Step 26d1: GU vs EU split of the NMPFams2 (metag) result.
# Queries: results/unk_filtered_aplakidou.faa (384,260 Aplakidou-filtered unknown 95% reps).
# Class per rep: results/fourway_classes_n88.tsv (class4). Hits: results/nmpfams2_metag_rep_hits.tsv
# (mmseqs search E <= 1e-3, >=30% id, -c 0.7 --cov-mode 0, vs family representatives).
# Read-only on every _n88 file.
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results"
OUT = R + "/nmpfams2_gu_eu.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")

q = [l[1:].split()[0] for l in open(R + "/unk_filtered_aplakidou.faa") if l[0] == ">"]
cls = {}
for i, l in enumerate(open(R + "/fourway_classes_n88.tsv")):
    if i == 0: continue
    f = l.split("\t", 2); cls[f[0]] = f[1]
hit = {l.split("\t", 1)[0] for l in open(R + "/nmpfams2_metag_rep_hits.tsv")}
miss = [k for k in q if k not in cls]
print("queries:", len(q), "| without a class:", len(miss))
if miss: raise SystemExit("ABORT: queries without a class")

n = collections.Counter(cls[k] for k in q)
h = collections.Counter(cls[k] for k in q if k in hit)
rows = [("class", "n_queries", "pct_of_queries", "represented", "not_represented", "pct_represented")]
for c in sorted(n) + ["total"]:
    a = len(q) if c == "total" else n[c]
    b = sum(h.values()) if c == "total" else h[c]
    rows.append((c, a, f"{100.0*a/len(q):.2f}", b, a - b, f"{100.0*b/a:.2f}"))
with open(OUT, "w") as fh:
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
for r in rows: print("\t".join(map(str, r)))
print("wrote", OUT)
