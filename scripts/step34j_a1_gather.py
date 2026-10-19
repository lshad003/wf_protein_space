# Step 34j A1: gene neighbourhoods (+-5 genes) of the strict members of the 4,509 recurrent EU families and of the
#   strict members of the 4,509 known length-matched controls, with neighbour annotation via 95% reps.
# Members: EU = results/step34h/strict_members.tsv (strict == 1) of every 95% rep of a recurrent family (rep to family
#   from results/step34c_strict_eu_nmpfams2.tsv). Controls = results/step34j/ctl_members.tsv with the Step 34h
#   eligibility rule (scripts/step34h_prep.py eligible(): rep-labelled contig Tiara-prokaryotic, fly90 == 0, LCA v2
#   not Eukaryota; computed contig additionally not fly >= 50% cover at >= 90% id and not frog >= 50% at >= 95%;
#   PAFs from results/step34g/paf or results/step34j/ctl_paf).
# Gene order, coordinates and strand: Prodigal headers of the per-assembly protein files (as step18a: catalog/
#   by_assembly_new, else shared Proteins/by_assembly); gene number = Prodigal order on the contig (contig_genenumber).
# Neighbour annotation (Step 29d method, n88 mapping): neighbour gene -> 95% rep (results/clusters_95_n88.tsv) and
#   50% cluster (results/clusters_50_n88.tsv); Pfam = results/pfam_hs/*.tblout (hmmsearch --cut_ga, Pfam-A, rep is the
#   target); eggNOG = results/eggnog2/*.emapper.annotations (KEGG_ko, KEGG_Module, Preferred_name, Description, COG).
# Orientation: offsets and strands are relative to the member (offset sign flipped when the member is on '-');
#   distance_bp = gap between member and neighbour CDS (negative = overlap).
# Outputs (refuse to overwrite): results/step34j/context_members.tsv, results/step34j/context_neighbours.tsv,
#   results/step34j/rep_annot.tsv (annotation of every neighbour 95% rep).
import os, gzip, glob, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34j/"
OLD = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"; NEW = W + "/catalog/by_assembly_new"
OM, ON, OA = D + "context_members.tsv", D + "context_neighbours.tsv", D + "rep_annot.tsv"
for p in (OM, ON, OA):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def contig(g):
    s, rest = g.split("__", 1); c, n = rest.rsplit("_", 1); return s, c, int(n)
fams = set()
with open(R + "step34c2c_strict_v2_families.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ir, ix = h.index("recurrent"), h.index("excluded")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ir] == "1" and f[ix] == "0": fams.add(f[0])
r2f = {}
with open(R + "step34c_strict_eu_nmpfams2.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[1] in fams: r2f[f[0]] = f[1]
print("recurrent families %d, their 95%% reps %d" % (len(fams), len(r2f)), flush=True)
if len(fams) != 4509 or set(r2f.values()) != fams: raise SystemExit("ABORT: family/rep mapping")
M = []   # (set, family, rep95, member)
n_eu_all = 0
with open(R + "step34h/strict_members.tsv") as fh:
    next(fh)
    for l in fh:
        r, m, u, s = l.rstrip("\n").split("\t")
        if r in r2f:
            n_eu_all += 1
            if s == "1": M.append(("EU", r2f[r], r, m))
print("EU members of recurrent families %d, strict %d" % (n_eu_all, len(M)), flush=True)
def paf_calls(p):
    iv, mt, bl, ql = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
    for l in open(p):
        f = l.split("\t"); k = f[0]
        ql[k] = int(f[1]); iv[k].append((int(f[2]), int(f[3]))); mt[k] += int(f[9]); bl[k] += int(f[10])
    out = {}
    for k, v in iv.items():
        v.sort(); c, cs, ce = 0, -1, -1
        for a, b in v:
            if a > ce: c += max(0, ce - cs); cs, ce = a, b
            else: ce = max(ce, b)
        out[k] = ((c + ce - cs) / ql[k], mt[k] / bl[k])
    return out
CM = [l.rstrip("\n").split("\t") for l in open(D + "ctl_members.tsv")]; H = CM[0]; CM = [dict(zip(H, r)) for r in CM[1:]]
fly, frog = {}, {}
for src, dd in (("step34g", R + "step34g/paf/"), ("step34j", D + "ctl_paf/")):
    stems = {m["stem"] for m in CM if m["paf_dir"] == src}
    for s in stems:
        for k, dct in (("fly", fly), ("frog", frog)):
            p = dd + "%s.%s.paf" % (s, k)
            if not os.path.exists(p): raise SystemExit("ABORT: missing " + p)
            for c, v in paf_calls(p).items(): dct[(s, c)] = v
def eligible(m):
    if m["label_source"] == "step34a_rep_label":
        return m["origin"] == "Tiara-prokaryotic" and m["fly90"] == "0" and m["lca_domain_v2"] != "Eukaryota"
    if m["origin"] != "Tiara-prokaryotic" or m["lca_domain_v2"] == "Eukaryota": return False
    k = (m["stem"], m["contig"]); a = fly.get(k, (0, 0)); b = frog.get(k, (0, 0))
    return not (a[0] >= 0.5 and a[1] >= 0.90) and not (b[0] >= 0.5 and b[1] >= 0.95)
nc = 0
for m in CM:
    if eligible(m): M.append(("control", m["rep"], m["rep"], m["member"])); nc += 1
print("control members %d, strict %d (%.2f%%)" % (len(CM), nc, 100.0 * nc / len(CM)), flush=True)
tc = {(contig(m)[0], contig(m)[1]) for _, _, _, m in M}
stems = sorted({s for s, _ in tc}); print("member contigs %d in %d stems" % (len(tc), len(stems)), flush=True)
genes = collections.defaultdict(dict)   # (stem, contig) -> n -> (start, end, strand, partial)
for s in stems:
    p = NEW + "/" + s + ".aa.fa.gz"
    if not os.path.exists(p): p = OLD + "/" + s + ".aa.fa.gz"
    with gzip.open(p, "rt") as fh:
        for l in fh:
            if l[0] != ">": continue
            f = l[1:].split(" # "); gid = f[0].strip(); c, n = gid.rsplit("_", 1)
            if (s, c) in tc:
                j = f[4].find("partial="); genes[(s, c)][int(n)] = (int(f[1]), int(f[2]), "+" if f[3] == "1" else "-", f[4][j + 8:j + 10])
print("contigs with gene coordinates %d of %d" % (len(genes), len(tc)), flush=True)
if len(genes) != len(tc): raise SystemExit("ABORT: contigs without headers")
gset = {"%s__%s_%d" % (s, c, n) for (s, c), d in genes.items() for n in d}
rep95, rep50 = {}, {}
for fn, dct in (("clusters_95_n88.tsv", rep95), ("clusters_50_n88.tsv", rep50)):
    with open(R + fn) as fh:
        for l in fh:
            r, m = l.rstrip("\n").split("\t")
            if m in gset: dct[m] = r
    print("%s: genes mapped %d of %d" % (fn, len(dct), len(gset)), flush=True)
reps = set(rep95.values())
pf = collections.defaultdict(set)
for f in glob.glob(R + "pfam_hs/*.tblout"):
    for l in open(f):
        if l[0] == "#": continue
        p = l.split()
        if p[0] in reps: pf[p[0]].add("%s|%s" % (p[3].split(".")[0], p[2]))
eg = {}
for f in glob.glob(R + "eggnog2/*.emapper.annotations"):
    for l in open(f):
        if l[0] == "#": continue
        p = l.rstrip("\n").split("\t")
        if p[0] in reps:
            eg[p[0]] = (p[11].replace("ko:", ""), p[13], p[8], p[7], p[6])
print("neighbour 95%% reps %d | with Pfam %d | with eggNOG %d" % (len(reps), len(pf), len(eg)), flush=True)
with open(OA, "w") as o:
    o.write("rep95\tpfam\tkegg_ko\tkegg_module\tpreferred_name\tdescription\tcog\n")
    for r in sorted(reps):
        e = eg.get(r, ("-",) * 5)
        o.write("%s\t%s\t%s\n" % (r, ",".join(sorted(pf.get(r, ()))) or "-", "\t".join(x if x else "-" for x in e)))
with open(OM, "w") as om, open(ON, "w") as on:
    om.write("set\tfamily\trep95\tmember\tstem\tcontig\tgene_n\tn_genes_contig\tstart\tend\tstrand\tpartial\tn_neighbours\n")
    on.write("set\tfamily\tmember\toffset\tneighbour\trel_strand\tdistance_bp\trep95\trep50\n")
    for st, fam, r, m in M:
        s, c, n = contig(m); g = genes[(s, c)]
        if n not in g: raise SystemExit("ABORT: member %s missing from its contig headers" % m)
        a, b, sd, pa = g[n]; k = 0
        for off in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5):
            if n + off not in g: continue
            a2, b2, sd2, _ = g[n + off]; ng = "%s__%s_%d" % (s, c, n + off)
            oo = off if sd == "+" else -off
            dist = (a2 - b - 1) if a2 > b else (a - b2 - 1)
            on.write("%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\n" % (st, fam, m, oo, ng, "same" if sd2 == sd else "opposite", dist, rep95.get(ng, "NA"), rep50.get(ng, "NA")))
            k += 1
        om.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n" % (st, fam, r, m, s, c, n, len(g), a, b, sd, pa, k))
print("wrote", OM, ON, OA)
