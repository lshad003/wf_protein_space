# Step 34c2b (report only; no set is changed): the 22 recurrent strict EU families whose majority genus is Drosophila
# (results/step34c2_recurrent_genus.tsv), and Eukaryota-labelled contigs across the strict set.
# Part 1, one row per (family, member contig): members via results/step34c_strict_eu_clusters30.tsv and
#   results/clusters_95_n88.tsv; Tiara class and probabilities from results/tiara/<stem>.tiara.tsv (contigs >= 1 kb);
#   contig LCA from <stem>_uniref50_lca.tsv[.gz] (rank, name, fragments, assigned, agreeing, support, lineage).
#   The LCA is per contig: gene_taxonomy_n88.tsv copies the contig label to every supported rep on it, so all genes on a
#   contig carry the same label; the supported reps on the contig (results/fourway_classes_n88.tsv) are listed with class.
#   dm6 alignments from the existing PAFs (results/step33g_paf and results/step34a_paf, <stem>.fly.paf; minimap2 -x asm10
#   --secondary=no -c). Only rep contigs were mapped (step33c_lists, step34a_lists); other member contigs are "not mapped".
#   Best alignment = PAF row with most matching bases: query coverage (qend-qstart)/qlen and identity col10/col11.
#   Fly rule as step34a: merged query coverage >= 0.50 AND sum(matches)/sum(block length) >= 0.95.
# Part 3: strict reps (step34a_rep_origin.tsv, Tiara-prokaryotic) whose contig is Eukaryota (a) by the domain column of
#   results/gene_taxonomy_n88.tsv and (b) by the lineage of the contig LCA, any token "<prefix>_Eukaryota" (step20c only
#   parsed "d_", and many stems write "-_Eukaryota", so (a) undercounts); by class and recurrent families touched.
# Output: results/step34c2b_drosophila_families.tsv; refuses to overwrite.
import os, re, gzip, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
FEC = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_scaffold_classify_mmseqs"
OUT = R + "step34c2b_drosophila_families.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
def contig(g):
    s, rest = g.split("__", 1); return (s, rest.rsplit("_", 1)[0])
fams = []
with open(R + "step34c2_recurrent_genus.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[0] == "recurrent_EU_family" and f[7] == "Drosophila": fams.append((f[1], f[8]))
print("recurrent families with majority genus Drosophila:", len(fams), flush=True)
if len(fams) != 22: raise SystemExit("ABORT: expected 22")
F = {f for f, _ in fams}
reps_of = collections.defaultdict(list); fam_of = {}
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t"); reps_of[a].append(b); fam_of[b] = a
rec = set()
with open(R + "step34c_strict_eu_families.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[8] == "1": rec.add(f[0])
need = {r for f in F for r in reps_of[f]}
mem = collections.defaultdict(list)
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if r in need: mem[r].append(m)
coh = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); coh[f[0]] = f[1]
orig = {}; cls = {}
with open(R + "step34a_rep_origin.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); orig[f[0]] = f[2]; cls[f[0]] = f[1]
mc = set()
for r in need:
    mc.add(contig(r))
    for m in mem[r]: mc.add(contig(m))
stems = sorted({s for s, _ in mc})
# supported reps on these contigs
oncontig = collections.defaultdict(list)
for r, c in cls.items():
    k = contig(r)
    if k in mc: oncontig[k].append(r)
# Tiara
tia = {}
for s in stems:
    p = R + "tiara/%s.tiara.tsv" % s
    if not os.path.exists(p): continue
    with open(p) as fh:
        next(fh)
        for l in fh:
            f = l.rstrip("\n").split("\t"); c = f[0].split()[0]
            if (s, c) in mc:
                tia[(s, c)] = (f[1], f[4], f[5], f[6], f[3], f[7], int(re.search(r"len=(\d+)", f[0]).group(1)))
# LCA, all contigs of all stems (part 1 detail + part 3 lineage domain)
def ldom(lin):
    for t in lin.split(";"):
        n = t.split("_", 1)[-1]
        if n in ("Bacteria", "Archaea", "Eukaryota", "Viruses"): return n
    return "unclassified"
lca = {}; ldm = {}; nolca = []
allstems = sorted({contig(r)[0] for r in cls})
for s in allstems:
    p = "%s/%s/%s_uniref50_lca.tsv" % (FEC, s, s)
    if not os.path.exists(p): p += ".gz"
    if not os.path.exists(p): nolca.append(s); continue
    for l in (gzip.open(p, "rt") if p.endswith(".gz") else open(p)):
        f = l.rstrip("\n").split("\t")
        if len(f) < 10: continue
        k = (s, f[0]); ldm[k] = ldom(f[9])
        if k in mc: lca[k] = (f[2], f[3], f[4], f[5], f[6], f[7], f[9])
print("stems in catalog: %d | without an LCA file: %d" % (len(allstems), len(nolca)), flush=True)
# fly PAFs
queried = set()
for s in stems:
    p = R + "step33c_lists/%s.tsv" % s
    if os.path.exists(p): queried |= {(s, l.split("\t")[0]) for l in open(p) if l.strip()}
    p = R + "step34a_lists/%s.txt" % s
    if os.path.exists(p): queried |= {(s, l.strip()) for l in open(p) if l.strip()}
aln = collections.defaultdict(list)
for s in stems:
    for d in ("step33g_paf", "step34a_paf"):
        p = R + "%s/%s.fly.paf" % (d, s)
        if not os.path.exists(p): continue
        for l in open(p):
            f = l.split("\t")
            if (s, f[0]) in mc: aln[(s, f[0])].append((int(f[1]), int(f[2]), int(f[3]), f[5], int(f[9]), int(f[10])))
def fly(k):
    if k not in queried: return ("not mapped (not a rep contig)", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
    v = aln.get(k)
    if not v: return ("no alignment", "NA", "NA", "NA", "NA", "NA", "NA", 0)
    ql = v[0][0]; b = max(v, key=lambda x: x[4])
    iv = sorted((x[1], x[2]) for x in v); c, cs, ce = 0, -1, -1
    for a, e in iv:
        if a > ce: c += max(0, ce - cs); cs, ce = a, e
        else: ce = max(ce, e)
    cov = (c + ce - cs) / ql; idn = sum(x[4] for x in v) / sum(x[5] for x in v)
    ok = cov >= 0.5 and idn >= 0.95
    why = "fly95" if ok else ("coverage < 50% and identity < 95%" if cov < 0.5 and idn < 0.95 else
                              ("coverage < 50%" if cov < 0.5 else "identity < 95%"))
    return (why, ql, b[3], "%.3f" % ((b[2] - b[1]) / ql), "%.3f" % (b[4] / b[5]), "%.3f" % cov, "%.3f" % idn, int(ok))
rows = []
for f, share in fams:
    for r in sorted(reps_of[f]):
        for m in sorted(set([r] + mem[r])):
            k = contig(m)
            t = tia.get(k); L = lca.get(k); fl = fly(k)
            ln = t[6] if t else (fl[1] if fl[1] != "NA" else "NA")
            rows.append([f, share, r, m, int(m == r), m.split("__")[0], coh.get(m.split("__")[0], "NA"), k[1], ln] +
                        (list(t[:6]) if t else ["<1 kb, not in Tiara"] + ["NA"] * 5) +
                        (list(L[:6]) + [ldm[k], L[6]] if L else ["no LCA row" if k[0] not in nolca else "no LCA file"] + ["NA"] * 7) +
                        [";".join("%s:%s" % (x.split("__")[1], cls[x]) for x in sorted(oncontig[k]))] + list(fl))
hdr = ["family_rep", "drosophila_share", "rep95", "member", "is_rep", "stem", "cohort", "contig", "contig_len", "tiara_class",
       "p_bac", "p_arc", "p_euk", "p_org", "p_unknown", "lca_rank", "lca_name", "lca_fragments", "lca_assigned", "lca_agree",
       "lca_support", "lca_domain_lineage", "lca_lineage", "supported_reps_on_contig", "fly_status", "fly_qlen",
       "fly_best_target", "fly_best_qcov", "fly_best_identity", "fly_rule_cov", "fly_rule_identity", "fly95"]
with open(OUT, "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for x in rows: fh.write("\t".join(map(str, x)) + "\n")
print("\nper family (rep contig details; member contig tallies)")
for f, share in fams:
    X = [x for x in rows if x[0] == f]
    for x in [x for x in X if x[4] == 1]:
        print("%s  Drosophila share %s | rep contig %s len %s | Tiara %s (bac %s euk %s) | LCA %s %s (frag %s, assigned %s, agree %s, support %s) | fly: %s best qcov %s id %s; rule cov %s id %s"
              % (f, share, x[7], x[8], x[9], x[10], x[12], x[15], x[16], x[17], x[18], x[19], x[20], x[24], x[27], x[28], x[29], x[30]))
        print("    supported reps on rep contig: %s" % x[23])
    mm = [x for x in X if x[4] == 0]
    print("    members %d (contigs %d); cohorts %s; Tiara %s; LCA name %s; fly %s" % (len(mm), len({(x[5], x[7]) for x in mm}),
          dict(collections.Counter(x[6] for x in mm)), dict(collections.Counter(x[9] for x in mm)),
          dict(collections.Counter(x[16] for x in mm).most_common(4)), dict(collections.Counter(x[24] for x in mm))))
print("\nsummary over the 22 rep contigs: fly status", dict(collections.Counter(x[24] for x in rows if x[4] == 1)))
print("rep contig LCA names", dict(collections.Counter(x[16] for x in rows if x[4] == 1)))
# Part 3
gtd = {}
with open(R + "gene_taxonomy_n88.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); gtd[f[0]] = f[4]
print("\nPart 3: strict reps on Eukaryota contigs (strict = Tiara-prokaryotic; %d stems have no LCA file, their reps count as no label)" % len(nolca))
print("%-5s %9s %22s %22s" % ("class", "strict", "gene_tax domain col", "lineage domain"))
cnt = collections.Counter(); a = collections.Counter(); b = collections.Counter()
eu_euk_a = set(); eu_euk_b = set()
for r, o in orig.items():
    if o != "Tiara-prokaryotic": continue
    c = cls[r]; cnt[c] += 1
    if gtd.get(r) == "Eukaryota": a[c] += 1; eu_euk_a.add(r)
    if ldm.get(contig(r)) == "Eukaryota": b[c] += 1; eu_euk_b.add(r)
for c in ("K", "KWP", "GU", "EU", "all"):
    n = sum(cnt.values()) if c == "all" else cnt[c]
    x = sum(a.values()) if c == "all" else a[c]; y = sum(b.values()) if c == "all" else b[c]
    print("%-5s %9d %12d (%6.3f%%) %12d (%6.3f%%)" % (c, n, x, 100.0 * x / n, y, 100.0 * y / n))
labelled = collections.Counter(cls[r] for r, o in orig.items() if o == "Tiara-prokaryotic" and contig(r) in ldm)
print("strict reps whose contig has an LCA row:", dict(labelled))
for lab, S in (("gene_tax domain col", eu_euk_a), ("lineage domain", eu_euk_b)):
    fr = {fam_of[r] for r in S if r in fam_of}
    print("recurrent families with >= 1 strict EU rep on a Eukaryota contig (%s): %d of %d (all families: %d)"
          % (lab, len(fr & rec), len(rec), len(fr)))
print("wrote", OUT, "rows:", len(rows))
