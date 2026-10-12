# Step 26d4: exact GU/EU union of NMPFams2 metag and iso representation (replaces the summed bound in 26d3).
# Queries: results/unk_filtered_aplakidou.faa (384,260 Aplakidou-filtered unknown 95% reps).
# Class per rep: results/fourway_classes_n88.tsv (class4). Hits: results/nmpfams2_metag_rep_hits.tsv and
# results/nmpfams2_iso_rep_hits.tsv (mmseqs search E <= 1e-3, >=30% id, -c 0.7 --cov-mode 0, vs family
# representatives). A query is represented in a set if it has any row in that set's hits file, as in step26d1.
# Read-only on every _n88 file.
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results"
OUT = R + "/nmpfams2_union_gu_eu.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")

q = [l[1:].split()[0] for l in open(R + "/unk_filtered_aplakidou.faa") if l[0] == ">"]
cls = {}
for i, l in enumerate(open(R + "/fourway_classes_n88.tsv")):
    if i == 0: continue
    f = l.split("\t", 2); cls[f[0]] = f[1]
m = {l.split("\t", 1)[0] for l in open(R + "/nmpfams2_metag_rep_hits.tsv")}
s = {l.split("\t", 1)[0] for l in open(R + "/nmpfams2_iso_rep_hits.tsv")}
miss = [k for k in q if k not in cls]
print("queries:", len(q), "| without a class:", len(miss))
if miss: raise SystemExit("ABORT: queries without a class")

rows = [("class", "n_queries", "metag", "iso", "both", "metag_only", "iso_only", "union", "neither",
         "pct_union", "pct_neither")]
for c in sorted(set(cls[k] for k in q)) + ["total"]:
    Q = set(q) if c == "total" else {k for k in q if cls[k] == c}
    a, b = Q & m, Q & s
    u = len(a | b)
    rows.append((c, len(Q), len(a), len(b), len(a & b), len(a - b), len(b - a), u, len(Q) - u,
                 f"{100.0*u/len(Q):.2f}", f"{100.0*(len(Q)-u)/len(Q):.2f}"))
with open(OUT, "w") as fh:
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
for r in rows: print("\t".join(map(str, r)))
print("wrote", OUT)
