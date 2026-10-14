# Step 30b: revision of step30a (gene neighbourhood of the 15,580 recurrent environmental-only 30% families).
# Changes from step30a:
# 1. Unit of conservation is the sequence variant: each distinct 95% rep among a family's members
#    (results/clusters_95_n88.tsv, the clustering the step26e1 members came from) is counted once, using its most
#    common +-3 neighbour annotation set over its member genes (mode over members with >= 1 neighbour gene; ties
#    broken by larger set, then sorted names). +-1 uses its own mode. Conservation = variants carrying the
#    annotation / all variants of the family. Families with one variant are "single-variant context": their
#    neighbours are reported (with member fractions) but they are never called conserved.
#    Genera: variant reps' contig taxonomy from results/gene_taxonomy_n88.tsv (supported reps only); genus = name
#    when rank is genus, first word of the name for ranks below genus ("Candidatus X" kept as two words), else none.
# 2. No COG letters: annotations are Pfam domains (results/pfam_hs/*.tblout, hmmsearch --cut_ga) and eggNOG
#    Preferred_name (results/eggnog2/*.emapper.annotations) only.
# 3. An annotation counts toward conservation only if it occurs within +-3 of a member of >= 5 different families
#    (counted over the 15,580 and the baseline families together).
# 4. Eukaryotic loci: distinct 95% reps of all neighbour genes of the family looked up in gene_taxonomy_n88.tsv;
#    euk_frac = Eukaryota / reps with a classified domain (not "unclassified"); flagged if euk_frac >= 0.5 with
#    >= 1 classified rep. euk_frac over all reps with a row is reported too.
# 5. Baseline: known reps = class4 K in results/fourway_classes_n88.tsv with pfam = 1 in
#    results/rep_classes_n88.tsv. A K 95% cluster is one variant by construction, so a K "family" is the 50% cluster
#    (results/clusters_50_n88.tsv, clustered independently from all proteins) containing the K rep; its genes are
#    the members and their 95% reps the variants. A pool of 60,000 K reps (seed 30) is mapped to 50% clusters
#    (deduplicated), then 2,000 families are drawn to match the 15,580 by n_variants bin (1, 2, 3, 4, 5-9, >= 10).
#    Members are restricted to the 88 samples; unlike step26e1 they are not filtered to complete >= 35 aa.
#    Same conservation statistic; enrichment background is pooled within each set.
# Conserved = n_variants >= 2, variants_pm3 >= 2, frac_pm3 >= 0.5 (also 0.75, 0.9), enrichment >= 5, annotation in
#   >= 5 families; the same rule for unknown and baseline families.
# Outputs (refuses to overwrite; step30a files are kept): results/step30b_family_summary.tsv,
#   results/step30b_family_annotations.tsv, results/step30b_single_variant_context.tsv,
#   results/step30b_category_map.tsv, results/step30b_baseline_comparison.tsv
import os, re, glob, random, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results/"
OUTS = [R + "step30b_" + x for x in ("family_summary.tsv", "family_annotations.tsv", "single_variant_context.tsv",
                                     "category_map.tsv", "baseline_comparison.tsv")]
for p in OUTS:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
OFF = (-3, -2, -1, 1, 2, 3)
POS = "UHM636.41030__k141_332084_5"
POOL, NBASE, MINFAM = 60000, 2000, 5
BINS = ((1, 1, "1"), (2, 2, "2"), (3, 3, "3"), (4, 4, "4"), (5, 9, "5-9"), (10, 10 ** 9, ">=10"))
binof = lambda n: next(b for lo, hi, b in BINS if lo <= n <= hi)
rng = random.Random(30)

# ---- unknown families (as step30a)
fams = [l[1:].split()[0] for l in open(R + "prostt5_main_eu15580.faa") if l[0] == ">"]
if len(fams) != 15580 or len(set(fams)) != 15580: raise SystemExit("ABORT: family set is not 15,580")
F = set(fams)
S88 = set()
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: S88.add(l.split("\t", 1)[0])
samp = lambda g: g.split("__", 1)[0]
animal = lambda g: samp(g).split(".", 1)[0]
mem = collections.defaultdict(list)
for l in open(R + "unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in F and samp(m) in S88: mem[r].append(m)
if len(mem) != 15580: raise SystemExit("ABORT: %d families have no retained member" % (15580 - len(mem)))
print("unknown families %d, members %d" % (len(mem), sum(map(len, mem.values()))), flush=True)

# ---- baseline pool: K reps with Pfam -> 50% clusters
K = set()
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        p = l.split("\t", 2)
        if p[1] == "K": K.add(p[0])
KP = set()
with open(R + "rep_classes_n88.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ip = h.index("pfam")
    for l in fh:
        p = l.rstrip("\n").split("\t")
        if p[0] in K and p[ip] == "1": KP.add(p[0])
print("K reps %d, with pfam = 1 %d" % (len(K), len(KP)), flush=True)
del K
pool = set(rng.sample(sorted(KP), POOL)); del KP
seed50 = {}
for l in open(R + "clusters_50_n88.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if m in pool: seed50[m] = r
R50 = collections.defaultdict(list)
for k, r in seed50.items(): R50[r].append(k)
print("pool %d K reps -> %d found in clusters_50_n88 -> %d distinct 50%% clusters" % (len(pool), len(seed50), len(R50)), flush=True)
cmem = collections.defaultdict(list)
for l in open(R + "clusters_50_n88.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in R50 and samp(m) in S88: cmem[r].append(m)
print("candidate baseline members %d" % sum(map(len, cmem.values())), flush=True)

def nb(g, off):  # as step29d
    stem, rest = g.split("__", 1); i = rest.rfind("_")
    try: n = int(rest[i + 1:])
    except ValueError: return None
    return "%s__%s_%d" % (stem, rest[:i], n + off)

# ---- member -> 95% rep (pass 1), then draw the matched baseline
need = {m for ms in mem.values() for m in ms} | {m for ms in cmem.values() for m in ms}
rep95 = {}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m in need: rep95[m] = r
print("members mapped to 95%% reps %d of %d" % (len(rep95), len(need)), flush=True)
nvar = lambda ms: len({rep95.get(m, m) for m in ms})
ubin = collections.Counter(binof(nvar(mem[q])) for q in fams)
cbin = collections.defaultdict(list)
for r in sorted(cmem): cbin[binof(nvar(cmem[r]))].append(r)
quota = {b: round(NBASE * ubin[b] / 15580) for _, _, b in BINS}
quota[max(quota, key=quota.get)] += NBASE - sum(quota.values())
base = []
print("n_variants bin: unknown families | candidate K 50% families | quota | drawn")
for _, _, b in BINS:
    got = rng.sample(cbin[b], min(quota[b], len(cbin[b]))); base += got
    print("  %-5s %6d %7d %5d %5d%s" % (b, ubin[b], len(cbin[b]), quota[b], len(got), "  SHORT" if len(got) < quota[b] else ""))
bmem = {r: cmem[r] for r in base}; del cmem
SETS = (("unknown", fams, mem), ("baseline", sorted(base), bmem))

# ---- neighbour genes -> 95% rep (pass 2)
need = {nb(m, o) for _, qs, M in SETS for q in qs for m in M[q] for o in OFF} - {None}
need = {n for n in need if n not in rep95}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m in need: rep95[m] = r
print("neighbour slots looked up %d" % len(need), flush=True)
nbrep = lambda g: rep95.get(g) if g else None

# ---- annotations (Pfam + eggNOG Preferred_name; no COG letters)
reps = set(rep95.values())
ann = collections.defaultdict(set)
for f in sorted(glob.glob(R + "pfam_hs/*.tblout")):
    for l in open(f):
        if l[0] == "#": continue
        p = l.split()
        if p[0] in reps: ann[p[0]].add("Pfam:" + p[2])
EGG = sorted(glob.glob(R + "eggnog2/*.emapper.annotations")); hdr = None
for f in EGG:
    for l in open(f):
        if l.startswith("#query"):
            h = l.lstrip("#").rstrip("\n").split("\t")
            if hdr is None: hdr = h; iP = h.index("Preferred_name"); print("eggNOG files %d, Preferred_name column %d" % (len(EGG), iP + 1))
            elif h != hdr: raise SystemExit("ABORT: eggNOG header differs in " + f)
            continue
        if l[0] == "#": continue
        p = l.rstrip("\n").split("\t")
        if p[0] in reps and p[iP] not in ("-", ""): ann[p[0]].add("eggNOG:" + p[iP])
print("95%% reps in play %d, with annotation %d" % (len(reps), sum(1 for r in reps if ann.get(r))), flush=True)

# ---- taxonomy of reps
BELOW = {"species", "strain", "species group", "species subgroup", "subspecies", "subgenus", "forma", "varietas",
         "isolate", "serotype", "serogroup", "biotype", "genotype", "forma specialis", "morph", "pathogroup"}
def genus(rank, name):
    if rank == "genus": return name
    if rank in BELOW:
        w = name.split()
        if not w or w[0] == "uncultured": return None
        return " ".join(w[:2]) if w[0] == "Candidatus" and len(w) > 1 else w[0]
    return None
tax = {}
with open(R + "gene_taxonomy_n88.tsv") as fh:
    next(fh)
    for l in fh:
        p = l.rstrip("\n").split("\t")
        if p[0] in reps: tax[p[0]] = (genus(p[2], p[3]), p[4])
print("reps with a taxonomy row %d" % len(tax), flush=True)

# ---- per family, per set
def mode(sets):
    c = collections.Counter(sets)
    return max(c, key=lambda s: (c[s], len(s), tuple(sorted(s)))) if c else frozenset()
FS = {}
for st, qs, M in SETS:
    for q in qs:
        byv = collections.defaultdict(list); anyann = set(); nbreps = set(); mw = 0
        for m in M[q]:
            s3, s1, has = set(), set(), False
            for o in OFF:
                r = nbrep(nb(m, o))
                if r is None: continue
                has = True; nbreps.add(r); d = ann.get(r, ()); s3.update(d)
                if abs(o) == 1: s1.update(d)
            mw += has; anyann |= s3
            byv[rep95.get(m, m)].append((has, frozenset(s3), frozenset(s1)))
        v3, v1, vnb = {}, {}, 0
        for v, xs in byv.items():
            w = [x for x in xs if x[0]]; vnb += bool(w)
            v3[v] = mode([x[1] for x in w]); v1[v] = mode([x[2] for x in w])
        mc3 = collections.Counter(a for xs in byv.values() for x in xs for a in x[1])
        FS[(st, q)] = dict(nm=len(M[q]), na=len({animal(m) for m in M[q]}), mw=mw, v3=v3, v1=v1, vnb=vnb, mc3=mc3,
                           anyann=anyann, nbreps=nbreps)
famocc = collections.Counter(a for s in FS.values() for a in s["anyann"])
bgv = {st: collections.Counter(a for (s2, q), s in FS.items() if s2 == st for v in s["v3"].values() for a in v) for st, _, _ in SETS}
nvtot = {st: sum(len(s["v3"]) for (s2, q), s in FS.items() if s2 == st) for st, _, _ in SETS}
print("variants: unknown %d, baseline %d | annotations in >= %d families: %d of %d" % (
    nvtot["unknown"], nvtot["baseline"], MINFAM, sum(1 for n in famocc.values() if n >= MINFAM), len(famocc)), flush=True)

# ---- category map (Pfam names and eggNOG preferred names only)
CMAP = [
 ("CRISPR/Cas and other defence", r"^Cas_|^Cas\d|CRISPR|^Cse\d|^Csm\d|^Cmr\d|^Csn\d|^Csy\d|^Csd\d|^Csx\d|^cas\d|^cas[A-E]$|^cs[a-z]\d"),
 ("CRISPR/Cas and other defence", r"^HsdR|^HsdM|^HsdS|Methylase_S|^RE_|Res_subunit|N6_N4_Mtase|DNA_methylase|^Abi|AbiEii|BREX|Gabija|Zorya|DISARM|^RelE|^ParE|^HicA|^HicB|^PIN(_|$)|^YafQ|^MazF|^PemK|^Txe"),
 ("phage/mobile", r"[Pp]hage|Terminase|Portal|Capsid|Baseplate|Holin|Prophage|Sipho|Tnp|Transposase|^DDE_|^rve|Integrase|Resolvase|Recombinase|Relaxase|^Mob[A-Z_]|^Tra[A-Z]|^Tail_|Phage_int|^int[A-Z]?$|^xer[CD]$"),
 ("transport/membrane", r"ABC_tran|ABC_membrane|BPD_transp|^MFS|[Tt]ransp(?!eptidase)|[Pp]ermease|Porin|^OEP|TonB|SBP_bac|Peripla_BP|Sugar_tr|^AA_permease|^Na_H|Channel|^MscS|^ECF|^Trk|^OmpA|^Plug$|^(mal|opp|dpp|liv|pst|pot|fec|fhu|feo|btu|mdt|acr|ton)[A-Z]$"),
 ("cell division/envelope", r"^Fts|^MinC|^MinD|^MinE|DivIVA|^ZapA|^SepF|^MreB|^RodA|Transpeptidase|^PBP|Peptidase_S11|^Mur_|^Mur[A-G]|^Lpx|^LptD|^Wza|Polysacc_synt|NLPC_P60|Amidase_[23]|^SpoVE|^LysM|^SLH|^(fts|min|mur|mre|rod|lpx|pbp)[A-Z]$"),
 ("transcription/regulation", r"^HTH|_HTH|GntR|LysR|TetR|MarR|AraC|LacI|Sigma70|^Sigma|Response_reg|^HisKA|Trans_reg_C|PadR|ArsR|^Crp$|^RHH|^Fur$|LuxR|Rrf2|AsnC|DeoR|^RNA_pol"),
 ("DNA replication/repair", r"^DNA_pol|^DnaB|^DnaG|^DnaA|^RecA|^Rec[FORJNQ]|^Uvr[ABCD]|^MutS|^MutL|DNA_ligase|^Ruv[ABC]|Topoisom|^Gyr|^SSB|^RNase_H|Toprim|^PolC|^(rec|uvr|mut|dna|ruv|gyr|top|ssb|lig)[A-Z]$"),
 ("metabolism", r"^adh_|Aldedh|Aldo_ket_red|Aminotran|Radical_SAM|Glyco_hydro|Pyr_redox|FAD_binding|NAD_binding|Oxidored|^Fer4|Rieske|^COX|Cytochrom|ATP-synt|TPP_enzyme|Beta_elim_lyase|^UbiA|Amidohydro|^HAD|Epimerase|Aldolase|Enolase|^PFK|Gp_dh|AMP-binding|Biotin|Carboxyl|Abhydrolase|^PALP|Ketoacyl|^GSDH|Dehydrogenase|_DH(_|$)"),
 ("secretion", r"^T[1-9]SS|^Sec[ADEFGY]|SecD_SecF|^YajC|^Tat[ABC]|^Gsp[A-Z]|Secretin|^VirB|^Hcp|^VgrG|T4SS|T6SS|T2SS|^(sec|tat|gsp)[A-Z]$"),
]
with open(OUTS[3], "w") as o:
    o.write("category\tregex\tnote\n")
    for c, rx in CMAP: o.write("%s\t%s\tcase-sensitive; Pfam name or eggNOG Preferred_name; COG letters not used\n" % (c, rx))
    o.write("other\t(no match, or matches of more than one category)\tambiguous matches are other\n")
CR = [(c, re.compile(rx)) for c, rx in CMAP]
def category(a):
    name = a.split(":", 1)[1]; hits = {c for c, rx in CR if rx.search(name)}
    return hits.pop() if len(hits) == 1 else ("other (ambiguous)" if hits else "other")

# ---- per family statistics and outputs
summ = {}
oa = open(OUTS[1], "w"); osv = open(OUTS[2], "w")
oa.write("set\tfamily_rep\tannotation\tn_variants\tvariants_pm3\tfrac_pm3\tvariants_pm1\tfrac_pm1\tn_members\tmembers_pm3\tbg_variant_frac\tenrichment\tn_families_with_annotation\tconserved50\tconserved75\tconserved90\tcategory\n")
osv.write("set\tfamily_rep\tn_members\tn_animals\tmembers_with_neighbour\tneighbour_euk_frac\teuk_flag\tcategories_in_context\tmodal_pm3_set(member_frac_pm3)\n")
for (st, q), s in FS.items():
    nv = len(s["v3"]); c3 = collections.Counter(a for v in s["v3"].values() for a in v); c1 = collections.Counter(a for v in s["v1"].values() for a in v)
    rows = []
    for x in sorted(c3):
        f3, f1 = c3[x] / nv, c1.get(x, 0) / nv; bg = bgv[st][x] / nvtot[st]; e = f3 / bg
        c = [nv >= 2 and c3[x] >= 2 and f3 >= t and e >= 5 and famocc[x] >= MINFAM for t in (0.5, 0.75, 0.9)]
        rows.append((x, f3, f1, e, c))
        oa.write("%s\t%s\t%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%d\t%.4g\t%.2f\t%d\t%s\t%s\t%s\t%s\n" % (st, q, x, nv, c3[x], f3, c1.get(x, 0), f1,
                 s["nm"], s["mc3"][x], bg, e, famocc[x], *["yes" if v else "no" for v in c], category(x)))
    cons = sorted([r for r in rows if r[4][0]], key=lambda r: (-r[1], -r[3], r[0]))
    mapped = [category(r[0]) for r in cons if not category(r[0]).startswith("other")]
    gen = [tax[v][0] for v in s["v3"] if v in tax and tax[v][0]]
    vdom = [tax[v][1] for v in s["v3"] if v in tax and tax[v][1] != "unclassified"]
    nt = [tax[r][1] for r in s["nbreps"] if r in tax]; ncl = [d for d in nt if d != "unclassified"]; ne = ncl.count("Eukaryota")
    ef = ne / len(ncl) if ncl else float("nan"); flag = bool(ncl) and ef >= 0.5
    sv = nv == 1
    summ[(st, q)] = dict(nm=s["nm"], na=s["na"], nv=nv, vnb=s["vnb"], mw=s["mw"], ngv=len(gen), ngen=len(set(gen)),
        gen=";".join(sorted(set(gen))[:10]) or "-", nnb=len(s["nbreps"]), nt=len(nt), ncl=len(ncl), ne=ne, ef=ef,
        efa=ne / len(nt) if nt else float("nan"), flag=flag, vef=vdom.count("Eukaryota") / len(vdom) if vdom else float("nan"),
        c50=len(cons), c75=sum(r[4][1] for r in rows), c90=sum(r[4][2] for r in rows),
        prim="single-variant context" if sv else (mapped[0] if mapped else ("other" if cons else "none")),
        cats=";".join(sorted(set(mapped))) or "-",
        top=";".join("%s(%.2f/%.2f/%.1f)" % (r[0], r[1], r[2], r[3]) for r in cons[:8]) or "-")
    if sv:
        (v,) = s["v3"]; ctx = sorted(s["v3"][v])
        ccat = sorted({category(x) for x in ctx} - {"other", "other (ambiguous)"})
        osv.write("%s\t%s\t%d\t%d\t%d\t%.3f\t%s\t%s\t%s\n" % (st, q, s["nm"], s["na"], s["mw"], ef, "yes" if flag else "no",
                  ";".join(ccat) or "-", ";".join("%s(%.2f)" % (x, s["mc3"][x] / s["nm"]) for x in ctx) or "-"))
        summ[(st, q)]["cats"] = ";".join(ccat) or "-"
oa.close(); osv.close()
with open(OUTS[0], "w") as o:
    o.write("set\tfamily_rep\tn_members\tn_animals\tn_variants\tvariant_class\tvariants_with_neighbour\tmembers_with_neighbour\tvariants_with_genus\tn_genera\tgenera(first10)\tneighbour_reps\tneighbour_reps_with_tax\tneighbour_reps_classified\tneighbour_reps_eukaryota\tneighbour_euk_frac\tneighbour_euk_frac_all\teuk_flag\tvariant_euk_frac\tn_conserved50\tn_conserved75\tn_conserved90\tprimary_category\tcategories\ttop_conserved(frac_pm3/frac_pm1/enrichment)\n")
    for (st, q), s in summ.items():
        o.write("%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%s\t%.3f\t%d\t%d\t%d\t%s\t%s\t%s\n" % (st, q, s["nm"], s["na"],
                s["nv"], "single-variant" if s["nv"] == 1 else "multi-variant", s["vnb"], s["mw"], s["ngv"], s["ngen"], s["gen"], s["nnb"], s["nt"],
                s["ncl"], s["ne"], s["ef"], s["efa"], "yes" if s["flag"] else "no", s["vef"], s["c50"], s["c75"], s["c90"], s["prim"], s["cats"], s["top"]))

# ---- report
U = {q: s for (st, q), s in summ.items() if st == "unknown"}; B = {q: s for (st, q), s in summ.items() if st == "baseline"}
print("\n== unknown families: single-variant %d, multi-variant %d (>= 3 variants %d)" % (
    sum(s["nv"] == 1 for s in U.values()), sum(s["nv"] >= 2 for s in U.values()), sum(s["nv"] >= 3 for s in U.values())))
print("   single-variant with >= 1 annotation in the modal pm3 set: %d" % sum(1 for q, s in U.items() if s["nv"] == 1 and any(FS[("unknown", q)]["v3"].values())))
print("   genera: families with >= 1 variant with genus %d; with >= 2 genera %d" % (sum(s["ngv"] > 0 for s in U.values()), sum(s["ngen"] >= 2 for s in U.values())))
print("   euk-flagged %d (of %d with a classified neighbour rep)" % (sum(s["flag"] for s in U.values()), sum(s["ncl"] > 0 for s in U.values())))
print("\n== conserved neighbourhood (>= 1 conserved annotation; variants_pm3 >= 2, enrichment >= 5, annotation in >= %d families)" % MINFAM)
print("   %-10s %-14s %8s %8s %8s %8s | %s" % ("set", "variants", "families", "c50", "c75", "c90", "c50 excluding euk-flagged (families)"))
for nm, D in (("unknown", U), ("baseline", B)):
    for k in (2, 3):
        X = [s for s in D.values() if s["nv"] >= k]; Y = [s for s in X if not s["flag"]]
        print("   %-10s >= %-11d %8d %8d %8d %8d | %d of %d" % (nm, k, len(X), sum(s["c50"] > 0 for s in X), sum(s["c75"] > 0 for s in X),
              sum(s["c90"] > 0 for s in X), sum(s["c50"] > 0 for s in Y), len(Y)))
with open(OUTS[4], "w") as o:
    o.write("set\tvariant_bin\tfamilies\tconserved50\tconserved50_frac\tconserved75\tconserved90\tfamilies_not_euk_flagged\tconserved50_not_euk_flagged\tconserved50_frac_not_euk_flagged\n")
    print("\n== baseline comparison by n_variants bin (conserved50 / families)")
    for lo, hi, b in BINS + ((2, 10 ** 9, ">=2"), (3, 10 ** 9, ">=3")):
        line = []
        for nm, D in (("unknown", U), ("baseline", B)):
            X = [s for s in D.values() if lo <= s["nv"] <= hi]; Y = [s for s in X if not s["flag"]]
            c, cy = sum(s["c50"] > 0 for s in X), sum(s["c50"] > 0 for s in Y)
            o.write("%s\t%s\t%d\t%d\t%.4f\t%d\t%d\t%d\t%d\t%.4f\n" % (nm, b, len(X), c, c / len(X) if X else float("nan"), sum(s["c75"] > 0 for s in X),
                    sum(s["c90"] > 0 for s in X), len(Y), cy, cy / len(Y) if Y else float("nan")))
            line.append("%s %d/%d (%.1f%%)" % (nm, c, len(X), 100 * c / len(X) if X else float("nan")))
        print("   %-5s %s" % (b, " | ".join(line)))
for lab, ex in (("all", False), ("excluding euk-flagged", True)):
    print("\n== primary category, unknown multi-variant families with a conserved annotation (%s)" % lab)
    for c, n in collections.Counter(s["prim"] for s in U.values() if s["nv"] >= 2 and s["c50"] and not (ex and s["flag"])).most_common(): print("   %-32s %d" % (c, n))
    print("   any category among conserved annotations:")
    for c, n in collections.Counter(c for s in U.values() if s["nv"] >= 2 and s["c50"] and not (ex and s["flag"]) for c in s["cats"].split(";") if c != "-").most_common(): print("   %-32s %d" % (c, n))
    print("   single-variant context, any category in the modal neighbour set (%s):" % lab)
    for c, n in collections.Counter(c for s in U.values() if s["nv"] == 1 and not (ex and s["flag"]) for c in s["cats"].split(";") if c != "-").most_common(): print("   %-32s %d" % (c, n))
print("\n== positive control %s" % POS)
s = U[POS]; f = FS[("unknown", POS)]
print("   members %d, animals %d, variants %d (%s), class %s, members with neighbour %d, genera %s, euk_frac %.3f, flag %s" % (
    s["nm"], s["na"], s["nv"], ";".join(f["v3"]), s["prim"], s["mw"], s["gen"], s["ef"], s["flag"]))
for v, st3 in f["v3"].items():
    print("   modal pm3 set: %s" % "; ".join("%s(%.2f, %s)" % (x, f["mc3"][x] / s["nm"], category(x)) for x in sorted(st3)))
    print("   modal pm1 set: %s" % "; ".join(sorted(f["v1"][v])))
print("done")
