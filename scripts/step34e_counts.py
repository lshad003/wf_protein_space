# Step 34e counts: subsets of results/count_matrix_primary_n88.tsv.gz (all 88 columns, unchanged values) for the gene
# sets analysed in Step 34e: strict_final (results/step34a_rep_origin_v2_final.tsv strict_final = 1; 1,359,273) and
# strict_final EU (41,114). Normalisation is done later within each set (step34e_hostfactor.R).
# Outputs (refuse to overwrite): results/step34e_counts_strict_final.tsv.gz, results/step34e_counts_strict_final_EU.tsv.gz
import os, gzip
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OA, OE = R + "step34e_counts_strict_final.tsv.gz", R + "step34e_counts_strict_final_EU.tsv.gz"
for p in (OA, OE):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
A, E = set(), set()
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, ifn = h.index("class4"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ifn] == "1":
            A.add(f[0])
            if f[ic] == "EU": E.add(f[0])
print("strict_final %d | EU %d" % (len(A), len(E)))
if len(A) != 1359273 or len(E) != 41114: raise SystemExit("ABORT: set sizes")
na = ne = 0
with gzip.open(R + "count_matrix_primary_n88.tsv.gz", "rt") as fi, gzip.open(OA, "wt") as oa, gzip.open(OE, "wt") as oe:
    hd = fi.readline(); oa.write(hd); oe.write(hd)
    for l in fi:
        r = l.split("\t", 1)[0]
        if r in A: oa.write(l); na += 1
        if r in E: oe.write(l); ne += 1
print("rows written: strict_final %d, EU %d" % (na, ne))
if na != len(A) or ne != len(E): raise SystemExit("ABORT: rows missing from count matrix")
