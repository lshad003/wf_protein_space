# Step 34r: per-representative x sample matrices from the identity-filtered recount (job 30125808,
#   scripts/step34r_remap_id95.sh). Rows = the 2,067,011 supported reps (results/supported_reps_95_n88.txt) in the row
#   order of results/count_matrix_primary_n88.tsv.gz; columns = its 88 samples in its order. One matrix per count column of
#   results/step34r/counts/<S>.counts.tsv.gz: col2 (primary mapped, as step16q col2), c_id95, c_id90, c_id95_q1.
# Checks: every sample file present; per-sample sum of col2 over all genes equals new_col2 in
#   results/step34r/summary/<S>.summary.tsv; c_id95 <= c_id90 <= col2 and c_id95_q1 <= c_id95 per cell.
# Outputs: results/step34r/matrix_{col2,c_id95,c_id90,c_id95_q1}_n88.tsv.gz, results/step34r/matrix_colsums_n88.tsv;
#   refuses to overwrite.
import os, gzip, sys
import numpy as np, pandas as pd
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; O = R + "step34r/"
NAMES = ["col2", "c_id95", "c_id90", "c_id95_q1"]
OUTS = [O + "matrix_%s_n88.tsv.gz" % n for n in NAMES] + [O + "matrix_colsums_n88.tsv"]
for p in OUTS:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
with gzip.open(R + "count_matrix_primary_n88.tsv.gz", "rt") as fh:
    samples = fh.readline().rstrip("\n").split("\t")[1:]
    genes = [l.split("\t", 1)[0] for l in fh]
if len(samples) != 88 or len(genes) != 2067011: raise SystemExit("ABORT: old matrix shape %d x %d" % (len(genes), len(samples)))
if set(genes) != set(open(R + "supported_reps_95_n88.txt").read().split()): raise SystemExit("ABORT: rows differ from supported list")
gi = pd.Index(genes)
M = np.zeros((4, len(genes), len(samples)), dtype=np.int64)
cs = []
for j, s in enumerate(samples):
    f = O + "counts/%s.counts.tsv.gz" % s
    if not os.path.exists(f): raise SystemExit("ABORT: missing " + f)
    d = pd.read_csv(f, sep="\t", header=None, names=["gene", "col2", "col3"] + NAMES[1:], dtype={"gene": str})
    sm = open(O + "summary/%s.summary.tsv" % s).read().rstrip("\n").split("\t")
    if int(d["col2"].sum()) != int(sm[7]): raise SystemExit("ABORT: %s col2 sum %d vs summary %s" % (s, d["col2"].sum(), sm[7]))
    if not ((d.c_id95 <= d.c_id90).all() and (d.c_id90 <= d.col2).all() and (d.c_id95_q1 <= d.c_id95).all()):
        raise SystemExit("ABORT: column order violated in " + s)
    k = gi.get_indexer(d["gene"]); ok = k >= 0
    for i, n in enumerate(NAMES): M[i, k[ok], j] = d[n].values[ok]
    cs.append([s, int(d["col2"].sum()), int((~ok).sum()), int(d.loc[~ok, "col2"].sum())] + [int(M[i, :, j].sum()) for i in range(4)])
    print("read", s, flush=True)
for i, n in enumerate(NAMES):
    df = pd.DataFrame(M[i], index=gi, columns=samples); df.index.name = "gene"
    df.to_csv(OUTS[i], sep="\t", compression={"method": "gzip", "compresslevel": 4})
    print("wrote", OUTS[i], "total", int(M[i].sum()), flush=True)
with open(OUTS[4], "w") as o:
    o.write("sample\tcol2_all_genes\tgenes_outside_supported\tcol2_outside_supported\t" + "\t".join(NAMES) + "\n")
    for r in cs: o.write("\t".join(map(str, r)) + "\n")
print("genes with col2 outside the supported rows, summed over samples: %d reads" % sum(r[3] for r in cs))
print("DONE step34r_matrices")
