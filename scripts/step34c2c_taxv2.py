# Step 34c2c part 1: domain parser fix for the contig taxonomy (step20c), and the summaries that read its domain column.
# step20c took the domain only from a "d_<name>" lineage token; stems whose lineage writes the domain as "-_<name>"
# (e.g. "-_cellular organisms;-_Bacteria") got "unclassified". v2 accepts any "<prefix>_<name>" token with name in
# Bacteria / Archaea / Eukaryota / Viruses. Same source, stems and rows as step20c_gene_tax_n88.py; rank and name
# unchanged (checked against the old file row by row).
# Output: results/gene_taxonomy_n88_v2.tsv (old file kept); results/step34c2c_summary.tsv (part 1 rows; part 2 appends).
# Scripts reading the domain column (grep of scripts/ for gene_taxonomy): step20c_gene_tax(_n88).py (writer; README
# Step 18), step30b_eu15580_neighbours.py (euk flag), step33e_tiara_vs_kraken.py (Tiara vs LCA), step27b_label_hits.py
# (contig domain), step34c2b_drosophila.py (part 3, already reported both). step34c2_genus.py reads rank/name only;
# step25c/25d only cite the file. Summaries are recomputed here from v2 with each script's own rule; no full reruns.
import os, re, gzip, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
FEC = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_scaffold_classify_mmseqs"
OUT = R + "gene_taxonomy_n88_v2.tsv"; SUM = R + "step34c2c_summary.tsv"
for p in (OUT, SUM):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
DOM = ("Bacteria", "Archaea", "Eukaryota", "Viruses")
def dom_old(lin):
    for t in lin.split(";"):
        if t.startswith("d_"): return t[2:]
    return "unclassified"
def dom_v2(lin):
    for t in lin.split(";"):
        n = t.split("_", 1)[-1]
        if n in DOM: return n
    return "unclassified"
summ = []
def S(sec, metric, old, new): summ.append((sec, metric, str(old), str(new))); print("  %-58s old %-22s new %s" % (metric, old, new), flush=True)
# ---- v2 file
cls = {}
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh: p = l.rstrip("\n").split("\t"); cls[p[0]] = p[1]
stems = (open(W + "/metadata/wf22_stems_44.txt").read().split() + open(W + "/metadata/wf23_stems.txt").read().split()
         + open(W + "/metadata/wf24_stems.txt").read().split())
def contig_of(g):
    s, rest = g.split("__", 1); return s + "__" + rest.rsplit("_", 1)[0]
need = collections.defaultdict(set)
for g in cls: need[g.split("__", 1)[0]].add(contig_of(g))
tax = {}; have = 0; tokfmt = collections.Counter()
for s in stems:
    f = "%s/%s/%s_uniref50_lca.tsv" % (FEC, s, s)
    if not os.path.exists(f): f += ".gz"
    if not os.path.exists(f): continue
    have += 1; fmt = collections.Counter()
    for l in (gzip.open(f, "rt") if f.endswith(".gz") else open(f)):
        p = l.rstrip("\n").split("\t")
        if len(p) < 10: continue
        k = s + "__" + p[0]
        if k in need[s]: tax[k] = (p[2], p[3], dom_old(p[9]), dom_v2(p[9]))
        m = re.search(r"(^|;)([^;_]*)_(Bacteria|Eukaryota|Archaea|Viruses)(;|$)", p[9])
        if m: fmt[m.group(2) + "_"] += 1
    tokfmt[",".join(sorted(fmt)) or "none"] += 1
print("stems with an LCA file: %d of %d | domain-token prefixes per stem: %s" % (have, len(stems), dict(tokfmt)), flush=True)
old = {}
with open(R + "gene_taxonomy_n88.tsv") as fh:
    next(fh)
    for l in fh: p = l.rstrip("\n").split("\t"); old[p[0]] = p
new = {}; bad = 0
with open(OUT, "w") as o:
    o.write("rep\tclass4\trank\tname\tdomain\n")
    for g, c in cls.items():
        t = tax.get(contig_of(g))
        if t is None: continue
        new[g] = (c, t[0], t[1], t[3])
        o.write("%s\t%s\t%s\t%s\t%s\n" % (g, c, t[0], t[1], t[3]))
        q = old.get(g)
        if q is None or q[1:5] != [c, t[0], t[1], t[2]]: bad += 1
print("v2 rows %d | old rows %d | rows whose rank/name/old-rule domain differ from the old file: %d" % (len(new), len(old), bad), flush=True)
if bad or set(new) != set(old): raise SystemExit("ABORT: v2 does not reproduce the old file under the old rule")
tr = collections.Counter((old[g][4], new[g][3]) for g in new)
print("domain transitions old -> v2:", {"%s->%s" % k: v for k, v in sorted(tr.items()) if k[0] != k[1]}, flush=True)
dO = {g: old[g][4] for g in new}; dN = {g: new[g][3] for g in new}
# ---- README Step 18 (step20c summary): domain by class among genes with contig taxonomy
print("\nREADME Step 18 (step20c): domain share by class, genes with contig taxonomy")
for c in ("K", "KWP", "GU", "EU"):
    G = [g for g in new if new[g][0] == c]
    for d in ("unclassified",) + DOM:
        a = sum(dO[g] == d for g in G); b = sum(dN[g] == d for g in G)
        S("README Step 18", "%s n=%d %s" % (c, len(G), d), "%d (%.1f%%)" % (a, 100.0 * a / len(G)), "%d (%.1f%%)" % (b, 100.0 * b / len(G)))
# ---- step33e part 2: Tiara (contigs >= 3 kb) x contig domain
GRP = {"bacteria": "prok", "archaea": "prok", "prokarya": "prok", "eukarya": "euk", "organelle": "euk", "unknown": "unknown"}
S88 = []
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: S88.append(l.split("\t", 1)[0])
tia = {}
for s in S88:
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 2)
            if int(re.search(r"len=(\d+)", f[0]).group(1)) >= 3000: tia[s + "__" + f[0].split()[0]] = GRP[f[1]]
print("\nstep33e: Tiara group x contig LCA domain, contigs >= 3 kb")
for lab, D in (("old", dO), ("v2", dN)):
    dom = {}
    for g, d in D.items():
        if d != "unclassified": dom[contig_of(g)] = d
    xt = collections.Counter((tia[k], dom[k]) for k in dom if k in tia)
    ag = xt[("prok", "Bacteria")] + xt[("prok", "Archaea")] + xt[("euk", "Eukaryota")]
    den = sum(xt[(g, d)] for g in ("prok", "euk") for d in ("Bacteria", "Archaea", "Eukaryota"))
    le = sum(xt[(g, "Eukaryota")] for g in ("prok", "euk", "unknown"))
    pe = sum(xt[("prok", d)] for d in ("Bacteria", "Archaea", "Eukaryota"))
    ee = sum(xt[("euk", d)] for d in ("Bacteria", "Archaea", "Eukaryota"))
    if lab == "old": o33 = (len(dom), sum(xt.values()), ag, den, xt, le, pe, ee)
    else: n33 = (len(dom), sum(xt.values()), ag, den, xt, le, pe, ee)
for nm, i in (("contigs with a domain (any length)", 0), ("contigs >= 3 kb with Tiara call and domain", 1)):
    S("step33e", nm, o33[i], n33[i])
for g in ("prok", "euk", "unknown"):
    for d in DOM: S("step33e", "Tiara %s x LCA %s" % (g, d), o33[4][(g, d)], n33[4][(g, d)])
f = lambda t: "%d/%d (%.2f%%)" % (t[2], t[3], 100.0 * t[2] / t[3])
S("step33e", "agreement prok-B/A + euk-E over prok/euk x B/A/E", f(o33), f(n33))
g = lambda t: "%d/%d (%.2f%%)" % (t[4][("euk", "Eukaryota")], t[5], 100.0 * t[4][("euk", "Eukaryota")] / t[5])
S("step33e", "LCA-Eukaryota contigs called euk by Tiara", g(o33), g(n33))
h = lambda t: "%d/%d (%.2f%%)" % (t[4][("prok", "Eukaryota")], t[6], 100.0 * t[4][("prok", "Eukaryota")] / t[6])
S("step33e", "Tiara-prok contigs with LCA Eukaryota (of B/A/E)", h(o33), h(n33))
# ---- step27b: Basidiobolus-vs-catalog targets, contig domain and bacterial flag (egg OR contig Bacteria)
print("\nstep27b: Basidiobolus hit targets")
rows = [l.rstrip("\n").split("\t") for l in open(R + "basid_vs_catalog_labelled.tsv")]
H = {c: i for i, c in enumerate(rows[0])}; rows = rows[1:]
tg = {}
for r in rows: tg[r[H["target"]]] = r[H["egg_lineage"]]
chk = sum(1 for r in rows if r[H["contig_domain"]] != dO.get(r[H["target"]], "NA"))
print("  labelled rows whose contig_domain differs from old gene_taxonomy:", chk)
if chk: raise SystemExit("ABORT: step27b labels do not reproduce from the old file")
for lab, D in (("old", dO), ("v2", dN)):
    td = collections.Counter(D.get(t, "NA") for t in tg)
    bac = {t for t, e in tg.items() if e == "Bacteria" or D.get(t) == "Bacteria"}
    co = sum(1 for t, e in tg.items() if e != "Bacteria" and D.get(t) == "Bacteria")
    bq = sum(1 for r in rows if r[H["is_best_hit"]] == "yes" and r[H["target"]] in bac)
    nq = sum(1 for r in rows if r[H["is_best_hit"]] == "yes")
    if lab == "old": o27 = (dict(td), len(bac), co, bq, nq)
    else: n27 = (dict(td), len(bac), co, bq, nq)
S("step27b", "targets by contig domain", ";".join("%s:%d" % kv for kv in sorted(o27[0].items())), ";".join("%s:%d" % kv for kv in sorted(n27[0].items())))
S("step27b", "targets bacterial (egg OR contig)", o27[1], n27[1])
S("step27b", "targets bacterial by contig only", o27[2], n27[2])
S("step27b", "queries whose best hit is bacterial", "%d (%.2f%%)" % (o27[3], 100.0 * o27[3] / o27[4]), "%d (%.2f%%)" % (n27[3], 100.0 * n27[3] / n27[4]))
# ---- step30b: euk flag of the 15,580 unknown families (neighbour reps +-3, euk_frac >= 0.5 over classified reps)
print("\nstep30b: euk flag")
OFF = (-3, -2, -1, 1, 2, 3)
fams = [l[1:].split()[0] for l in open(R + "prostt5_main_eu15580.faa") if l[0] == ">"]
F = set(fams); S8 = set(S88)
mem = collections.defaultdict(list)
for l in open(R + "unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in F and m.split("__", 1)[0] in S8: mem[r].append(m)
def nb(g, off):
    stem, rest = g.split("__", 1); i = rest.rfind("_")
    try: n = int(rest[i + 1:])
    except ValueError: return None
    return "%s__%s_%d" % (stem, rest[:i], n + off)
slots = {nb(m, o) for ms in mem.values() for m in ms for o in OFF} - {None}
rep95 = {}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m in slots: rep95[m] = r
nbr = {q: {rep95[x] for m in mem[q] for o in OFF for x in [nb(m, o)] if x in rep95} for q in fams}
fs = {}
with open(R + "step30b_family_summary.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        p = l.rstrip("\n").split("\t")
        if p[0] == "unknown": fs[p[1]] = dict(zip(h, p))
def flags(D):
    out = {}
    for q in fams:
        nt = [D[r] for r in nbr[q] if r in D]; ncl = [d for d in nt if d != "unclassified"]
        out[q] = (bool(ncl) and ncl.count("Eukaryota") / len(ncl) >= 0.5, bool(ncl), ncl.count("Eukaryota") / len(ncl) if ncl else float("nan"))
    return out
fo, fn = flags(dO), flags(dN)
mis = sum(1 for q in fams if ("yes" if fo[q][0] else "no") != fs[q]["euk_flag"])
print("  families %d | old-rule flags differing from step30b_family_summary.tsv: %d" % (len(fams), mis), flush=True)
if mis: raise SystemExit("ABORT: old euk flags do not reproduce step30b")
cnt = lambda fl: "%d (of %d with a classified neighbour rep)" % (sum(v[0] for v in fl.values()), sum(v[1] for v in fl.values()))
S("step30b", "euk-flagged families of 15,580", cnt(fo), cnt(fn))
for lab, sel in (("multi-variant (>= 2)", lambda q: int(fs[q]["n_variants"]) >= 2), (">= 3 variants", lambda q: int(fs[q]["n_variants"]) >= 3)):
    Q = [q for q in fams if sel(q)]; c50 = [q for q in Q if int(fs[q]["n_conserved50"]) > 0]
    e = lambda fl: "c50 %d of %d not flagged (flagged %d of %d families)" % (sum(1 for q in c50 if not fl[q][0]), sum(1 for q in Q if not fl[q][0]), sum(fl[q][0] for q in Q), len(Q))
    S("step30b", lab + ": conserved50 excluding euk-flagged", e(fo), e(fn))
sv = {}
with open(R + "step30b_single_variant_context.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        p = l.rstrip("\n").split("\t")
        if p[0] == "unknown": sv[p[1]] = dict(zip(h, p))
d84 = [q for q, v in sv.items() if "CRISPR/Cas and other defence" in v["categories_in_context"].split(";")]
S("step30b", "defence-context single-variant families (%d): euk-flagged" % len(d84), sum(fo[q][0] for q in d84), sum(fn[q][0] for q in d84))
P = "UHM636.41030__k141_332084_5"
S("step30b", "Cse2 family euk_frac / flag", "%.3f / %s" % (fo[P][2], fo[P][0]), "%.3f / %s" % (fn[P][2], fn[P][0]))
with open(SUM, "w") as o:
    o.write("section\tmetric\told\tnew\n")
    for r in summ: o.write("\t".join(r) + "\n")
print("wrote", OUT, "and", SUM, "(part 1)")
