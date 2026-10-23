# Step 34s check: direction and class mix of the c_id95 WF22 A_blocked STP1717.1 genes (padj < 0.05), top 10 by |logFC|.
# Inputs: results/step34s/c_id95/step34e_limma/strict_final_WF22_A_blocked_treatmentSTP1717_1.csv.gz,
#   results/step34a_rep_origin_v2_final.tsv (class4, origin, lca_domain_v2).
# Outputs in results/step34s/stp1717_check/: sig459_classes.tsv (summary), sig459_genes.tsv (all significant), top10_abs_logfc.tsv
import csv, gzip
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; O = W + "/results/step34s/stp1717_check/"
L = W + "/results/step34s/c_id95/step34e_limma/strict_final_WF22_A_blocked_treatmentSTP1717_1.csv.gz"
sig = []
with gzip.open(L, "rt") as fh:
    r = csv.reader(fh); next(r)
    for f in r:
        if float(f[5]) < 0.05: sig.append(f)
ids = {f[0] for f in sig}; cl = {}
with open(W + "/results/step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[0] in ids: cl[f[0]] = (f[1], f[2], f[7])
print("significant", len(sig), "classed", len(cl))
with open(O + "sig459_genes.tsv", "w") as o:
    o.write("rep\tclass4\torigin\tlca_domain_v2\tlogFC\tAveExpr\tP.Value\tadj.P.Val\n")
    for f in sorted(sig, key=lambda f: -abs(float(f[1]))):
        o.write("\t".join([f[0], *cl[f[0]], f[1], f[2], f[4], f[5]]) + "\n")
S = {}
for f in sig:
    k = cl[f[0]][0]; d = "up" if float(f[1]) > 0 else "down"
    S.setdefault(k, {"up": 0, "down": 0})[d] += 1
with open(O + "sig459_classes.tsv", "w") as o:
    o.write("class4\tup\tdown\ttotal\n")
    for k in ["K", "KWP", "GU", "EU"]:
        s = S.get(k, {"up": 0, "down": 0}); o.write("%s\t%d\t%d\t%d\n" % (k, s["up"], s["down"], s["up"] + s["down"]))
    o.write("all\t%d\t%d\t%d\n" % (sum(s["up"] for s in S.values()), sum(s["down"] for s in S.values()), len(sig)))
with open(O + "sig459_genes.tsv") as fh, open(O + "top10_abs_logfc.tsv", "w") as o:
    for i, l in enumerate(fh):
        if i <= 10: o.write(l)
