# Step 27c: contig context of Basidiobolus bacterial-candidate targets (step27b).
# Input: results/basid_bacterial_candidates.tsv (one best-hit row per query), filtered here to qcov >= 0.8 AND
#   tcov >= 0.8.
# Contig and gene position from the gene ID, as in Step 21 (step21b_neighbours_n88.py): stem__contig_n, contig =
#   everything before the last "_". Contig genes = every member in results/clusters_95_n88.tsv (col 2) on the same
#   stem__contig; each is mapped to its 95% representative (col 1). UHM586.41010 genes are absent from the n88
#   table, as in Step 21 n88, so its contigs have no neighbours here.
# Neighbour lineage: eggNOG of the neighbour's 95% rep, same rule as step27b (results/eggnog2, eggNOG_OGs levels):
#   Bacteria, Fungi, other Eukaryota, Archaea, root_only, conflict, none.
# Context per target, over annotated neighbours (lineage not "none"): "bacterial context" if > 50% are Bacteria;
#   "fungal context" if > 50% are Fungi or other Eukaryota; "mixed" otherwise; "no annotated neighbours" if none
#   (including single-gene contigs).
# Read-only on every _n88 file.
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results"
OUT = R + "/basid_candidates_context.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")

def contig(g):
    stem, rest = g.split("__", 1)
    i = rest.rfind("_")
    return stem + "__" + rest[:i], int(rest[i + 1:])

L = open(R + "/basid_bacterial_candidates.tsv").readline().rstrip("\n").split("\t")
C = {c: i for i, c in enumerate(L)}
cand = [l.rstrip("\n").split("\t") for i, l in enumerate(open(R + "/basid_bacterial_candidates.tsv")) if i > 0]
keep = [r for r in cand if float(r[C["qcov"]]) >= 0.8 and float(r[C["tcov"]]) >= 0.8]
print("candidate queries:", len(cand), "| qcov >= 0.8 and tcov >= 0.8:", len(keep))
byt = collections.defaultdict(list)
for r in keep: byt[r[C["target"]]].append(r)
ctg = {contig(t)[0] for t in byt}
print("targets:", len(byt), "| contigs:", len(ctg))

genes = collections.defaultdict(dict)  # contig -> {gene: rep}
for l in open(R + "/clusters_95_n88.tsv"):
    r, m = l.rstrip("\n").split("\t")[:2]
    c = m.rsplit("_", 1)[0]
    if c in ctg: genes[c][m] = r
need = {r for c in genes for r in genes[c].values()}
print("contig genes found:", sum(len(v) for v in genes.values()), "| distinct reps:", len(need))

egg = {}
nfiles = 0
for n in range(620):
    p = R + "/eggnog2/LsPS.%d.emapper.annotations" % n
    if not os.path.exists(p): continue
    nfiles += 1
    for l in open(p):
        if l[0] == "#": continue
        k = l.split("\t", 1)[0]
        if k not in need: continue
        f = l.rstrip("\n").split("\t")
        lv = {x.split("@", 1)[1].split("|", 1)[0] for x in f[4].split(",") if "@" in x}
        b, e = "2" in lv, ("2759" in lv or "4751" in lv)
        if b and e: lin = "conflict"
        elif b: lin = "Bacteria"
        elif "4751" in lv: lin = "Fungi"
        elif "2759" in lv: lin = "other Eukaryota"
        elif "2157" in lv: lin = "Archaea"
        else: lin = "root_only"
        egg[k] = (lin, f[5])
if nfiles != 619: raise SystemExit("ABORT: expected 619 eggNOG files")

HDR = ["target", "stem", "cohort", "treatment", "group", "supported", "target_egg_lineage", "target_max_annot_lvl",
       "n_queries", "strains", "queries", "min_pident", "min_qcov", "contig", "n_contig_genes", "n_neighbours",
       "n_annotated", "n_bacteria", "n_fungi_euk", "n_other", "context", "neighbours"]
rows = []
for t, rs in byt.items():
    c, pos = contig(t)
    nb = sorted(((g, rp) for g, rp in genes.get(c, {}).items() if g != t), key=lambda x: contig(x[0])[1])
    lab = [(g, egg.get(rp, ("none", "NA"))[0]) for g, rp in nb]
    a = [x for _, x in lab if x != "none"]
    nbac = sum(x == "Bacteria" for x in a)
    nfe = sum(x in ("Fungi", "other Eukaryota") for x in a)
    if not a: ctx = "no annotated neighbours"
    elif nbac > len(a) / 2.0: ctx = "bacterial context"
    elif nfe > len(a) / 2.0: ctx = "fungal context"
    else: ctx = "mixed"
    r0 = rs[0]
    rows.append([t, r0[C["stem"]], r0[C["cohort"]], r0[C["treatment"]], r0[C["group"]], r0[C["supported"]],
                 r0[C["egg_lineage"]], r0[C["max_annot_lvl"]], str(len(rs)),
                 ",".join(sorted({r[C["strain"]] for r in rs})), ",".join(sorted(r[C["query"]] for r in rs)),
                 "%.3f" % min(float(r[C["pident"]]) for r in rs), "%.3f" % min(float(r[C["qcov"]]) for r in rs),
                 c, str(len(genes.get(c, {}))), str(len(nb)), str(len(a)), str(nbac), str(nfe),
                 str(len(a) - nbac - nfe), ctx,
                 ",".join("%s:%s" % (g.rsplit("_", 1)[1], x) for g, x in lab) or "NA"])
order = ["bacterial context", "fungal context", "mixed", "no annotated neighbours"]
rows.sort(key=lambda r: (order.index(r[20]), -int(r[8]), r[0]))
with open(OUT, "w") as fh:
    fh.write("\t".join(HDR) + "\n")
    for r in rows: fh.write("\t".join(r) + "\n")

H = {c: i for i, c in enumerate(HDR)}
print("\ncontext            targets  queries  supported_targets")
for o in order:
    a = [r for r in rows if r[H["context"]] == o]
    print("%-24s %4d %8d %6d" % (o, len(a), sum(int(r[H["n_queries"]]) for r in a),
                                 sum(r[H["supported"]] == "yes" for r in a)))
print("targets on contigs with no other gene:", sum(r[H["n_neighbours"]] == "0" for r in rows))
print("targets with no contig genes in clusters_95_n88 (UHM586 or absent):", sum(r[H["n_contig_genes"]] == "0" for r in rows))
print("\nbacterial-context targets:")
S = ["target", "stem", "cohort", "treatment", "supported", "n_queries", "strains", "min_pident", "n_annotated",
     "n_bacteria", "target_max_annot_lvl"]
print("\t".join(S))
for r in rows:
    if r[H["context"]] == "bacterial context": print("\t".join(r[H[c]] for c in S))
print("\nbacterial-context by group (targets / queries):")
for g in ("inoculated", "control", "excluded"):
    a = [r for r in rows if r[H["context"]] == "bacterial context" and r[H["group"]] == g]
    print("  %-10s %d / %d" % (g, len(a), sum(int(r[H["n_queries"]]) for r in a)))
print("\nwrote", OUT)
