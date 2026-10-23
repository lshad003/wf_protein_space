# Step 36: genomic context (neighbouring genes) of the final recurrent EU families, with a K baseline.
# Replaces step30b (old 15,580-family set, not citable) and sits beside the Step 34j Part A context (+-5 genes,
#   50%-cluster Jaccard), which is not changed.
# Families: the 4,509 recurrent strict_final EU families (results/step34c2c_strict_v2_families.tsv, recurrent == 1,
#   excluded == 0; 95% reps from results/step34c_strict_eu_nmpfams2.tsv). Flags: member_validated = every 95% cluster
#   class_any == EU (results/step34h_member_classes.tsv, as step34j_c1); structure_confident =
#   structure_confident_annotated == 1 in results/step34j_structures_v2.tsv (ProtNLM-only names unannotated).
# Members: proteins of the families' 95% clusters on strict-eligible contigs (results/step34h/strict_members.tsv,
#   strict == 1), restricted to the 88 metagenomes of metadata/samples_88_treatment.tsv (UHM586.41010 dropped).
# K baseline: the 4,509 length-matched K controls of Step 34j (results/step34j/esmfold_input_labels.tsv, set control);
#   members = results/step34j/ctl_members.tsv with the Step 34j A1 eligibility rule (copied from step34j_a1_gather.py),
#   same 88 metagenomes.
# Three runs, one procedure: EU_family30 (all members of the 30% family; the EU table), EU_rep95 (members of the 95%
#   cluster of each family's Step 34j query, results/step34d_query_map.tsv) and K_rep95 (members of each control's own
#   95% cluster). The EU-vs-K comparison is EU_rep95 vs K_rep95.
# Window: genes within 2 positions of the member (offsets -2..+2 by Prodigal gene number) on the same contig and the
#   same strand. Gene number, coordinates and strand from the original Prodigal headers (catalog/by_assembly_new, else
#   the shared Proteins/by_assembly, as step9 and step34j_a1); files are streamed and only member contigs kept.
#   Offsets are oriented to the member (sign flipped when the member is on '-').
# Labels: neighbour gene -> 95% rep (results/clusters_95_n88.tsv) -> class: class4 of results/fourway_classes_n88.tsv
#   (supported reps; K, KWP, GU, EU), else class of results/rep_classes_n88.tsv (K, KWP, U). K: each Pfam-A domain of
#   the rep (results/pfam_hs/*.tblout, hmmsearch --cut_ga, Pfam 38.2) whose name does not start DUF or UPF ("Pfam:<name>");
#   a K rep with only DUF/UPF domains is DUF-only. KWP: eggNOG Preferred_name (results/eggnog2/*.emapper.annotations;
#   "eggNOG:<name>"). Uninformative: GU, EU, U, DUF-only, KWP with no Preferred_name, gene without a 95% rep.
# Conserved label (criterion): carried by the neighbours of >= 50% of the unit's member contigs (denominator = all
#   member contigs, with or without neighbours), in >= 3 contigs from >= 2 animals. Per conserved label: fraction of
#   contigs, modal oriented position (and its share of contigs), and the number of distinct 95% clusters of the unit
#   whose member contigs carry it. A unit's conserved context is "several_variants" if any conserved label is supported
#   by >= 2 distinct 95% clusters, else "one_variant" (strain-level locus). For the rep95 runs every unit is one 95%
#   cluster, so the split applies to EU_family30 only.
# Positive control: the Cse2 family UHM45.35758__k141_100597_236 must have a conserved cas Pfam label (Step 34f type
#   I-E names); otherwise the script writes no output and exits 2.
# Outputs (refuse to overwrite): results/step36/context_by_family.tsv, results/step36/conserved_labels.tsv,
#   results/step36/structure_vs_context.tsv, results/step36/member_windows.tsv (one row per member and neighbour).
import os, gzip, glob, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step36/"
OLD = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"; NEW = W + "/catalog/by_assembly_new"
OF, OL, OS, OW = D + "context_by_family.tsv", D + "conserved_labels.tsv", D + "structure_vs_context.tsv", D + "member_windows.tsv"
for p in (OF, OL, OS, OW):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
CSE2 = "UHM45.35758__k141_100597_236"
CAS = {"CRISPR_Cse1", "CRISPR_Cse2", "Cas_CT1975", "Cas_Cas5d", "CRISPR_assoc", "Cas_Cas1", "CRISPR_Cas2", "Cas_Cas2CT1978",
       "Cas3_HD", "DEAD", "Helicase_C"}
FRAC, MINC, MINA, WIN = 0.5, 3, 2, 2
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
def contig(g):
    s, rest = g.split("__", 1); c, n = rest.rsplit("_", 1); return s, c, int(n)
animal = lambda s: s.split(".", 1)[0]
S88 = {r["stem"] for r in tsv(W + "/metadata/samples_88_treatment.tsv")}
if len(S88) != 88 or any(s.startswith("UHM586") for s in S88): raise SystemExit("ABORT: sample set")
# ---- EU families, flags, reps
fams = {r["family_rep"] for r in tsv(R + "step34c2c_strict_v2_families.tsv") if r["recurrent"] == "1" and r["excluded"] == "0"}
r2f = {r["rep"]: r["family_rep"] for r in tsv(R + "step34c_strict_eu_nmpfams2.tsv") if r["family_rep"] in fams}
if len(fams) != 4509 or set(r2f.values()) != fams: raise SystemExit("ABORT: family/rep mapping")
mc = collections.defaultdict(list)
for r in tsv(R + "step34h_member_classes.tsv"):
    if r["family_rep"] in fams: mc[r["family_rep"]].append(r["class_any"])
mval = {f: int(bool(mc[f]) and all(c == "EU" for c in mc[f])) for f in fams}
st = {r["family_rep"]: r for r in tsv(R + "step34j_structures_v2.tsv") if r["set"] == "EU"}
conf = {f for f, r in st.items() if r["structure_confident_annotated"] == "1"}
qm = {r["family_rep"]: r["query"] for r in tsv(R + "step34d_query_map.tsv") if r["family_rep"] in fams}
print("EU families %d | 95%% reps %d | member-validated %d | structure-confident %d | queries %d"
      % (len(fams), len(r2f), sum(mval.values()), len(conf), len(qm)), flush=True)
if sum(mval.values()) != 4250 or len(conf) != 94 or len(qm) != 4509: raise SystemExit("ABORT: flag counts")
# ---- members: (run, unit, rep95, member)
UNITS = collections.defaultdict(list)   # (run, unit) -> [(rep95, member)]
n_all = n_strict = n_drop = 0
qrep = {q: f for f, q in qm.items()}
for r in tsv(R + "step34h/strict_members.tsv"):
    if r["rep"] not in r2f: continue
    n_all += 1
    if r["strict"] != "1": continue
    n_strict += 1
    if contig(r["member"])[0] not in S88: n_drop += 1; continue
    UNITS[("EU_family30", r2f[r["rep"]])].append((r["rep"], r["member"]))
    if r["rep"] in qrep: UNITS[("EU_rep95", qrep[r["rep"]])].append((r["rep"], r["member"]))
print("EU members of recurrent families %d | strict %d | dropped (not in the 88) %d" % (n_all, n_strict, n_drop), flush=True)
miss = [q for q in qm.values() if ("EU_rep95", qrep[q]) not in UNITS]
print("EU queries whose own 95%% cluster has no strict member in the 88: %d" % len(miss), flush=True)
# controls (eligibility as step34j_a1_gather.py)
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
ctl = {r["id"] for r in tsv(R + "step34j/esmfold_input_labels.tsv") if r["set"] == "control"}
CM = list(tsv(R + "step34j/ctl_members.tsv"))
fly, frog = {}, {}
for src, dd in (("step34g", R + "step34g/paf/"), ("step34j", R + "step34j/ctl_paf/")):
    for s in {m["stem"] for m in CM if m["paf_dir"] == src}:
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
nc = nd = 0
for m in CM:
    if m["rep"] not in ctl: raise SystemExit("ABORT: control member of an unknown rep " + m["rep"])
    if not eligible(m): continue
    nc += 1
    if m["stem"] not in S88: nd += 1; continue
    UNITS[("K_rep95", m["rep"])].append((m["rep"], m["member"]))
print("controls %d | members %d | eligible %d | dropped (not in the 88) %d | controls with an eligible member %d"
      % (len(ctl), len(CM), nc, nd, sum(1 for k in UNITS if k[0] == "K_rep95")), flush=True)
if len(ctl) != 4509: raise SystemExit("ABORT: control count")
# ---- Prodigal headers of member contigs
allm = {m for v in UNITS.values() for _, m in v}
tc = {contig(m)[:2] for m in allm}
genes = collections.defaultdict(dict)
for s in sorted({s for s, _ in tc}):
    p = NEW + "/" + s + ".aa.fa.gz"
    if not os.path.exists(p): p = OLD + "/" + s + ".aa.fa.gz"
    with gzip.open(p, "rt") as fh:
        for l in fh:
            if l[0] != ">": continue
            f = l[1:].split(" # "); c, n = f[0].strip().rsplit("_", 1)
            if (s, c) in tc: genes[(s, c)][int(n)] = (int(f[1]), int(f[2]), "+" if f[3] == "1" else "-")
print("member genes %d | member contigs %d | contigs with headers %d" % (len(allm), len(tc), len(genes)), flush=True)
if len(genes) != len(tc): raise SystemExit("ABORT: contigs without Prodigal headers")
win = {}   # member -> [(oriented offset, neighbour gene)]
for m in allm:
    s, c, n = contig(m); g = genes[(s, c)]
    if n not in g: raise SystemExit("ABORT: member not in its contig headers " + m)
    sd = g[n][2]; v = []
    for off in range(-WIN, WIN + 1):
        if off and n + off in g and g[n + off][2] == sd:
            v.append((off if sd == "+" else -off, "%s__%s_%d" % (s, c, n + off)))
    win[m] = v
need = {x for v in win.values() for _, x in v}
print("members with >= 1 same-strand neighbour within 2 genes: %d of %d | neighbour genes %d"
      % (sum(1 for v in win.values() if v), len(win), len(need)), flush=True)
# ---- neighbour labels
rep95 = {}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m in need: rep95[m] = r
reps = set(rep95.values())
print("neighbour genes with a 95%% rep %d of %d | distinct reps %d" % (len(rep95), len(need), len(reps)), flush=True)
cls = {}
for r in tsv(R + "fourway_classes_n88.tsv"):
    if r["rep"] in reps: cls[r["rep"]] = r["class4"]
n4 = len(cls)
for r in tsv(R + "rep_classes_n88.tsv"):
    if r["rep"] in reps and r["rep"] not in cls: cls[r["rep"]] = r["class"]
print("rep class from fourway_classes_n88 %d, from rep_classes_n88 %d, none %d" % (n4, len(cls) - n4, len(reps - set(cls))), flush=True)
pf = collections.defaultdict(set)
for f in glob.glob(R + "pfam_hs/*.tblout"):
    for l in open(f):
        if l[0] == "#": continue
        p = l.split()
        if p[0] in reps: pf[p[0]].add(p[2])
eg = {}
for f in glob.glob(R + "eggnog2/*.emapper.annotations"):
    for l in open(f):
        if l[0] == "#": continue
        p = l.rstrip("\n").split("\t")
        if p[0] in reps: eg[p[0]] = p[8]
def labels(gene):
    r = rep95.get(gene)
    if r is None: return "no_rep", []
    c = cls.get(r, "NA")
    if c == "K":
        x = sorted("Pfam:" + d for d in pf.get(r, ()) if not (d.startswith("DUF") or d.startswith("UPF")))
        return ("K", x) if x else ("DUF_only", [])
    if c == "KWP":
        n = eg.get(r, "-")
        return ("KWP", ["eggNOG:" + n]) if n not in ("-", "") else ("KWP_no_name", [])
    return c, []
LAB = {g: labels(g) for g in need}
print("neighbour genes by label class:", dict(collections.Counter(v[0] for v in LAB.values())), flush=True)
with open(OW + ".tmp", "w") as o:
    o.write("run\tunit\trep95\tmember\toffset\tneighbour\tneighbour_rep95\tneighbour_class\tlabels\n")
    for (run, u), v in sorted(UNITS.items()):
        for r, m in v:
            for off, x in win[m]:
                o.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (run, u, r, m, off, x, rep95.get(x, "NA"), LAB[x][0], ",".join(LAB[x][1]) or "-"))
# ---- conservation per unit
def score(v):
    ctg = collections.defaultdict(lambda: dict(reps=set(), lab=collections.defaultdict(set)))
    for r, m in v:
        s, c, _ = contig(m); k = (s, c); ctg[k]["reps"].add(r)
        for off, x in win[m]:
            for L in LAB[x][1]: ctg[k]["lab"][L].add(off)
    nC = len(ctg); out = []
    allL = collections.Counter(L for d in ctg.values() for L in d["lab"])
    for L, n in allL.items():
        ks = [k for k, d in ctg.items() if L in d["lab"]]
        an = {animal(k[0]) for k in ks}
        pos = collections.Counter(min(ctg[k]["lab"][L], key=abs) for k in ks)
        p, pn = max(pos.items(), key=lambda t: (t[1], -abs(t[0]), t[0]))
        nr = len({r for k in ks for r in ctg[k]["reps"]})
        ok = n / nC >= FRAC and n >= MINC and len(an) >= MINA
        out.append(dict(label=L, n_contigs=n, frac=n / nC, n_animals=len(an), modal_offset=p, modal_offset_frac=pn / n, n_rep95=nr, conserved=ok))
    out.sort(key=lambda d: (-d["frac"], -d["n_rep95"], d["label"]))
    return nC, ctg, out
rows, crow = [], []
for (run, u), v in sorted(UNITS.items()):
    nC, ctg, out = score(v)
    cons = [d for d in out if d["conserved"]]
    fam = u if run != "EU_rep95" else qrep[u]
    if run == "K_rep95": fam = "-"
    vc = "none" if not cons else ("several_variants" if any(d["n_rep95"] >= 2 for d in cons) else "one_variant")
    an = {animal(k[0]) for k in ctg}
    rows.append(dict(run=run, unit=u, family_rep=fam, member_validated=mval.get(fam, "NA") if run != "K_rep95" else "NA",
        structure_confident=int(fam in conf) if run != "K_rep95" else "NA", n_rep95=len({r for r, _ in v}), n_members=len(v),
        n_members_with_neighbour=sum(1 for _, m in v if win[m]), n_member_contigs=nC, n_animals=len(an),
        n_contigs_informative=sum(1 for d in ctg.values() if d["lab"]), callable=int(nC >= MINC and len(an) >= MINA),
        n_conserved_labels=len(cons), context_class=vc,
        conserved_labels=";".join("%s:%.2f:%+d:%d" % (d["label"], d["frac"], d["modal_offset"], d["n_rep95"]) for d in cons) or "-",
        top_labels=";".join("%s:%.2f" % (d["label"], d["frac"]) for d in out[:3]) or "-"))
    for d in cons: crow.append(dict(run=run, unit=u, family_rep=fam, n_member_contigs=nC, n_rep95_unit=len({r for r, _ in v}), **d))
# empty units (no strict member in the 88) still get a row
have = {(r["run"], r["unit"]) for r in rows}
for run, us in (("EU_family30", sorted(fams)), ("EU_rep95", sorted(qm.values())), ("K_rep95", sorted(ctl))):
    for u in us:
        if (run, u) not in have:
            fam = u if run == "EU_family30" else (qrep[u] if run == "EU_rep95" else "-")
            rows.append(dict(run=run, unit=u, family_rep=fam, member_validated=mval.get(fam, "NA") if run != "K_rep95" else "NA",
                structure_confident=int(fam in conf) if run != "K_rep95" else "NA", n_rep95=0, n_members=0, n_members_with_neighbour=0,
                n_member_contigs=0, n_animals=0, n_contigs_informative=0, callable=0, n_conserved_labels=0, context_class="none",
                conserved_labels="-", top_labels="-"))
# ---- Cse2 positive control
c = [d for d in crow if d["run"] == "EU_family30" and d["unit"] == CSE2]
cr = [r for r in rows if r["run"] == "EU_family30" and r["unit"] == CSE2][0]
print("\nCse2 family %s: member contigs %d, animals %d, 95%% clusters %d" % (CSE2, cr["n_member_contigs"], cr["n_animals"], cr["n_rep95"]))
for d in c: print("   conserved %-22s %d/%d contigs (%.2f), %d animals, modal offset %+d (%.2f), 95%% clusters %d"
                  % (d["label"], d["n_contigs"], cr["n_member_contigs"], d["frac"], d["n_animals"], d["modal_offset"], d["modal_offset_frac"], d["n_rep95"]))
casc = [d["label"] for d in c if d["label"].split(":", 1)[1] in CAS]
if not casc:
    os.remove(OW + ".tmp"); print("POSITIVE CONTROL FAILED: no conserved cas label; no outputs written"); raise SystemExit(2)
print("positive control passed: conserved cas labels", casc)
# ---- write
os.rename(OW + ".tmp", OW)
order = {"EU_family30": 0, "EU_rep95": 1, "K_rep95": 2}
rows.sort(key=lambda r: (order[r["run"]], r["unit"]))
H = list(rows[0].keys())
with open(OF, "w") as o:
    o.write("\t".join(H) + "\n")
    for r in rows: o.write("\t".join(str(r[k]) for k in H) + "\n")
HL = ["run", "unit", "family_rep", "n_member_contigs", "n_rep95_unit", "label", "n_contigs", "frac", "n_animals", "modal_offset",
      "modal_offset_frac", "n_rep95"]
with open(OL, "w") as o:
    o.write("\t".join(HL) + "\tsupport\n")
    for d in sorted(crow, key=lambda d: (order[d["run"]], d["unit"], -d["frac"], d["label"])):
        o.write("\t".join(("%.3f" % d[k]) if k in ("frac", "modal_offset_frac") else str(d[k]) for k in HL)
                + "\t%s\n" % ("several_variants" if d["n_rep95"] >= 2 else "one_variant"))
fr = {r["unit"]: r for r in rows if r["run"] == "EU_family30"}
with open(OS, "w") as o:
    o.write("family_rep\tquery\tmean_ca_plddt\thit_class\tconfident_db\tconfident_target\tconfident_target_name\tconfident_qtm\t"
            "n_rep95\tn_member_contigs\tn_animals\tcontext_class\tconserved_labels\ttop_labels\n")
    for f in sorted(conf, key=lambda f: (-float(st[f]["confident_qtm"]), f)):
        s, r = st[f], fr[f]
        o.write("\t".join(str(x) for x in (f, qm[f], s["mean_ca_plddt"], s["hit_class"], s["confident_db"], s["confident_target"],
                s["confident_target_name"].replace("\t", " "), s["confident_qtm"], r["n_rep95"], r["n_member_contigs"], r["n_animals"],
                r["context_class"], r["conserved_labels"], r["top_labels"])) + "\n")
# ---- summary
print("\nconserved context (>= %.0f%% of member contigs, >= %d contigs, >= %d animals; same strand, within %d genes)" % (FRAC * 100, MINC, MINA, WIN))
def summ(name, R_):
    n = len(R_); cb = sum(r["callable"] for r in R_); k = collections.Counter(r["context_class"] for r in R_)
    a = n - k["none"]
    print("  %-34s units %5d | callable (>= 3 contigs, >= 2 animals) %5d | any conserved %5d (%.2f%% of units, %.2f%% of callable)"
          " | one variant %d, several variants %d" % (name, n, cb, a, 100.0 * a / n, 100.0 * a / cb if cb else 0, k["one_variant"], k["several_variants"]))
E30 = [r for r in rows if r["run"] == "EU_family30"]
summ("EU_family30, all", E30)
summ("EU_family30, member-validated", [r for r in E30 if r["member_validated"] == 1])
summ("EU_family30, structure-confident", [r for r in E30 if r["structure_confident"] == 1])
summ("EU_rep95 (query's own 95% cluster)", [r for r in rows if r["run"] == "EU_rep95"])
summ("K_rep95 (control's own 95% cluster)", [r for r in rows if r["run"] == "K_rep95"])
print("wrote", OF, OL, OS, OW)
