# Step 34s check: per-sample reads of the 459 c_id95 WF22 STP1717.1 genes (sig459_genes.tsv) in the c_id95 strict_final
#   counts, WF22 columns only. Output: results/step34s/stp1717_check/sig459_per_sample.tsv (stem, animal, treatment,
#   egg_mass, month, run, total reads over the 459, number of the 459 with >= 1 read, share of WF22 library in the 459).
import gzip
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; O = W + "/results/step34s/stp1717_check/"
ids = {l.split("\t")[0] for l in open(O + "sig459_genes.tsv")} - {"rep"}
meta = [l.rstrip("\n").split("\t") for l in open(W + "/metadata/wf22_design.tsv")][1:]
run = {f[0]: f[2] for f in (l.rstrip("\n").split("\t") for l in open(W + "/results/batch_table_88.tsv"))}
with gzip.open(W + "/results/step34s/c_id95/step34e_counts_strict_final.tsv.gz", "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {s: h.index(s) for s in (m[0] for m in meta)}
    tot = dict.fromkeys(ix, 0); sub = dict.fromkeys(ix, 0); nz = dict.fromkeys(ix, 0)
    for l in fh:
        f = l.rstrip("\n").split("\t"); g = f[0] in ids
        for s, i in ix.items():
            c = int(f[i]); tot[s] += c
            if g: sub[s] += c; nz[s] += c > 0
with open(O + "sig459_per_sample.tsv", "w") as o:
    o.write("stem\tanimal\ttreatment\tegg_mass\tmonth\trun\treads_459\tgenes_ge1\tshare_pct\n")
    for m in meta:
        s = m[0]; o.write("\t".join(m + [run[s], str(sub[s]), str(nz[s]), "%.4f" % (100 * sub[s] / tot[s])]) + "\n")
