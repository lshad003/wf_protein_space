# Step 30a: gene neighbourhood of the 15,580 recurrent environmental-only 30% families (Step 28 query set),
# extending step29d (same +-3 offsets and neighbour-ID parsing) from 4 families to all 15,580.
# Families: headers of results/prostt5_main_eu15580.faa (abort unless 15,580); cross-checked, not re-derived,
#   against results/unk_direct30_eu_recurrence.tsv (n_animals >= 10 and n_cohorts >= 2).
# Members: results/unk_direct30_clusters.tsv (rep, member); members whose sample is not in
#   metadata/samples_88_treatment.tsv are dropped (UHM586.41010 excluded). Animal = sample prefix before ".".
# Neighbours: genes at offsets -3..-1, +1..+3 on the member's contig, mapped to 95% reps via results/clusters_95.tsv
#   (streamed; needed IDs only). Annotations of each neighbour rep:
#   Pfam: every domain in results/pfam_hs/*.tblout (hmmsearch --cut_ga), as "Pfam:<name>";
#   eggNOG: results/eggnog2/*.emapper.annotations (scripts/step8f_eggnog_v2.sh), Preferred_name as
#   "eggNOG:<name>" and each COG_category letter as "COG:<letter>"; no other eggNOG column is used.
# Per family and annotation: fraction of distinct animals (primary) and of members with it within +-3 and +-1.
# Background (all 15,580 families pooled):
#   bg_slot  = neighbour slots carrying the annotation / all neighbour slots with a gene (as specified);
#   bg_win3  = members whose +-3 window carries it / all members (same unit as the family fraction).
#   enrichment = animal_frac_pm3 / bg_slot; enrichment_win = animal_frac_pm3 / bg_win3.
# Conserved (as specified): animal_frac_pm3 >= 0.5 and enrichment >= 5; also at 0.75 and 0.9.
#   conserved_win: same with enrichment_win >= 5, reported alongside.
# Categories: regex map in results/eu15580_neighbour_category_map.tsv, applied to conserved annotations only;
#   an annotation matching more than one category is "other (ambiguous)". No function is assigned to the
#   family proteins themselves.
# HIT families (best E over PDB and CATH50 <= 1e-3, results/prostt5_main/eu15580_vs_{PDB,CATH50}.tsv): top hit
#   (lowest E, then highest bits) with PDB title from the Foldseek PDB_h header file, beside conserved annotations.
# Outputs (refuses to overwrite): results/eu15580_neighbour_summary.tsv, results/eu15580_neighbour_annotations.tsv,
#   results/eu15580_neighbour_category_map.tsv, results/eu15580_hit_vs_context.tsv
import os, re, glob, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results/"
FDB = "/bigdata/stajichlab/lshad003/software/foldseek_dbs"
OUTS = [R + x for x in ("eu15580_neighbour_summary.tsv", "eu15580_neighbour_annotations.tsv",
                         "eu15580_neighbour_category_map.tsv", "eu15580_hit_vs_context.tsv")]
for p in OUTS:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
OFF = (-3, -2, -1, 1, 2, 3)
POS = "UHM636.41030__k141_332084_5"

# ---- family set
fams = [l[1:].split()[0] for l in open(R + "prostt5_main_eu15580.faa") if l[0] == ">"]
if len(fams) != 15580 or len(set(fams)) != 15580: raise SystemExit("ABORT: family set is %d (%d unique), not 15,580" % (len(fams), len(set(fams))))
F = set(fams)
rec = set()
with open(R + "unk_direct30_eu_recurrence.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ia, ic = h.index("n_animals"), h.index("n_cohorts")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if int(f[ia]) >= 10 and int(f[ic]) >= 2: rec.add(f[0])
print("families: faa %d | recurrence rule %d | faa only %d | recurrence only %d" % (len(F), len(rec), len(F - rec), len(rec - F)))
for x in sorted(F - rec)[:20]: print("  faa only:", x)
for x in sorted(rec - F)[:20]: print("  recurrence only:", x)

# ---- members
S88 = set()
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: S88.add(l.split("\t", 1)[0])
mem = collections.defaultdict(list); drop = collections.Counter()
for l in open(R + "unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r not in F: continue
    s = m.split("__", 1)[0]
    if s not in S88: drop[s] += 1; continue
    mem[r].append(m)
print("members kept %d in %d families | dropped %d (%s)" % (sum(map(len, mem.values())), len(mem), sum(drop.values()), dict(drop)))
if len(mem) != 15580: raise SystemExit("ABORT: %d families have no retained member" % (15580 - len(mem)))
animal = lambda m: m.split("__", 1)[0].split(".", 1)[0]

def nb(g, off):  # as step29d
    stem, rest = g.split("__", 1); i = rest.rfind("_")
    try: n = int(rest[i + 1:])
    except ValueError: return None
    return "%s__%s_%d" % (stem, rest[:i], n + off)
need = {nb(m, o) for ms in mem.values() for m in ms for o in OFF} - {None}
print("neighbour slots needed %d" % len(need), flush=True)

# ---- neighbour -> 95% rep
rep = {}
with open(R + "clusters_95.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")[:2]
        if m in need: rep[m] = r
reps = set(rep.values())
print("neighbour genes found %d | distinct 95%% reps %d" % (len(rep), len(reps)), flush=True)

# ---- annotations
ann = collections.defaultdict(set)
for f in sorted(glob.glob(R + "pfam_hs/*.tblout")):
    for l in open(f):
        if l[0] == "#": continue
        p = l.split()
        if p[0] in reps: ann[p[0]].add("Pfam:" + p[2])
print("reps with Pfam %d" % sum(1 for r in reps if any(a.startswith("Pfam:") for a in ann.get(r, ()))), flush=True)
EGG = sorted(glob.glob(R + "eggnog2/*.emapper.annotations"))
print("eggNOG annotation files: %s/eggnog2/LsPS.<N>.emapper.annotations, %d files" % (R.rstrip("/"), len(EGG)))
hdr = None
for f in EGG:
    for l in open(f):
        if l.startswith("#query"):
            h = l.lstrip("#").rstrip("\n").split("\t")
            if hdr is None:
                hdr = h; print("eggNOG header:", " | ".join(h))
                iP, iC = h.index("Preferred_name"), h.index("COG_category")
                print("using columns %d Preferred_name, %d COG_category (1-based)" % (iP + 1, iC + 1), flush=True)
            elif h != hdr: raise SystemExit("ABORT: eggNOG header differs in " + f)
            continue
        if l[0] == "#": continue
        p = l.rstrip("\n").split("\t")
        if p[0] not in reps: continue
        if p[iP] not in ("-", ""): ann[p[0]].add("eggNOG:" + p[iP])
        for c in p[iC]:
            if c not in "-": ann[p[0]].add("COG:" + c)
print("reps with any annotation %d of %d" % (sum(1 for r in reps if ann.get(r)), len(reps)), flush=True)

# ---- per member windows and background
slots = 0; bg_slot = collections.Counter(); bg_w3 = collections.Counter(); bg_w1 = collections.Counter(); nmem = 0
fam_stats = {}
for q in fams:
    a3m, a1m, a3a, a1a = collections.Counter(), collections.Counter(), collections.defaultdict(set), collections.defaultdict(set)
    an_all, m_nb, an_nb = set(), 0, set()
    for m in mem[q]:
        A = animal(m); an_all.add(A); nmem += 1
        s3, s1, has = set(), set(), False
        for o in OFF:
            n = nb(m, o)
            if n not in rep: continue
            has = True; slots += 1; d = ann.get(rep[n], set())
            bg_slot.update(d); s3 |= d
            if abs(o) == 1: s1 |= d
        if has: m_nb += 1; an_nb.add(A)
        bg_w3.update(s3); bg_w1.update(s1)
        a3m.update(s3); a1m.update(s1)
        for x in s3: a3a[x].add(A)
        for x in s1: a1a[x].add(A)
    fam_stats[q] = (len(mem[q]), len(an_all), m_nb, len(an_nb), a3m, a1m, {k: len(v) for k, v in a3a.items()}, {k: len(v) for k, v in a1a.items()})
print("background: %d neighbour slots with a gene, %d members, %d distinct annotations" % (slots, nmem, len(bg_slot)), flush=True)

# ---- category map
CMAP = [
 ("CRISPR/Cas and other defence", "any", r"^Cas_|^Cas\d|CRISPR|^Cse\d|^Csm\d|^Cmr\d|^Csn\d|^Csy\d|^Csd\d|^Csx\d|^cas\d|^cs[a-z]\d"),
 ("CRISPR/Cas and other defence", "any", r"^HsdR|^HsdM|^HsdS|Methylase_S|^RE_|Res_subunit|N6_N4_Mtase|DNA_methylase|^Abi|AbiEii|BREX|Gabija|Zorya|DISARM|^RelE|^ParE|^HicA|^HicB|^PIN(_|$)|^YafQ|^MazF|^PemK|^Txe"),
 ("CRISPR/Cas and other defence", "COG", r"^V$"),
 ("phage/mobile", "any", r"[Pp]hage|Terminase|Portal|Capsid|Baseplate|Holin|Prophage|Sipho|Tnp|Transposase|^DDE_|^rve|Integrase|Resolvase|Recombinase|Relaxase|^Mob[A-Z_]|^Tra[A-Z]|^Tail_|Phage_int|^int[A-Z]?$|^xer[CD]$"),
 ("phage/mobile", "COG", r"^X$"),
 ("transport/membrane", "any", r"ABC_tran|ABC_membrane|BPD_transp|^MFS|[Tt]ransp(?!eptidase)|[Pp]ermease|Porin|^OEP|TonB|SBP_bac|Peripla_BP|Sugar_tr|^AA_permease|^Na_H|Channel|^MscS|^ECF|^Trk|^OmpA|^Plug$|^(mal|opp|dpp|liv|pst|pot|fec|fhu|feo|btu|mdt|acr|ton)[A-Z]$"),
 ("cell division/envelope", "any", r"^Fts|^MinC|^MinD|^MinE|DivIVA|^ZapA|^SepF|^MreB|^RodA|Transpeptidase|^PBP|Peptidase_S11|^Mur_|^Mur[A-G]|^Lpx|^LptD|^Wza|Polysacc_synt|NLPC_P60|Amidase_[23]|^SpoVE|^LysM|^SLH|^(fts|min|mur|mre|rod|lpx|pbp)[A-Z]$"),
 ("cell division/envelope", "COG", r"^[DM]$"),
 ("transcription/regulation", "any", r"^HTH|_HTH|GntR|LysR|TetR|MarR|AraC|LacI|Sigma70|^Sigma|Response_reg|^HisKA|Trans_reg_C|PadR|ArsR|^Crp$|^RHH|^Fur$|LuxR|Rrf2|AsnC|DeoR|^RNA_pol"),
 ("transcription/regulation", "COG", r"^K$"),
 ("DNA replication/repair", "any", r"^DNA_pol|^DnaB|^DnaG|^DnaA|^RecA|^Rec[FORJNQ]|^Uvr[ABCD]|^MutS|^MutL|DNA_ligase|^Ruv[ABC]|Topoisom|^Gyr|^SSB|^RNase_H|Toprim|^PolC|^(rec|uvr|mut|dna|ruv|gyr|top|ssb|lig)[A-Z]$"),
 ("DNA replication/repair", "COG", r"^L$"),
 ("metabolism", "any", r"^adh_|Aldedh|Aldo_ket_red|Aminotran|Radical_SAM|Glyco_hydro|Pyr_redox|FAD_binding|NAD_binding|Oxidored|^Fer4|Rieske|^COX|Cytochrom|ATP-synt|TPP_enzyme|Beta_elim_lyase|^UbiA|Amidohydro|^HAD|Epimerase|Aldolase|Enolase|^PFK|Gp_dh|AMP-binding|Biotin|Carboxyl|Abhydrolase|^PALP|Ketoacyl|^GSDH|Dehydrogenase|_DH(_|$)"),
 ("metabolism", "COG", r"^C$"),
 ("secretion", "any", r"^T[1-9]SS|^Sec[ADEFGY]|SecD_SecF|^YajC|^Tat[ABC]|^Gsp[A-Z]|Secretin|^VirB|^Hcp|^VgrG|T4SS|T6SS|T2SS|^(sec|tat|gsp)[A-Z]$"),
 ("secretion", "COG", r"^U$"),
]
with open(OUTS[2], "w") as o:
    o.write("category\tsource\tregex\tnote\n")
    for c, s, rx in CMAP: o.write("%s\t%s\t%s\t%s\n" % (c, s, rx, "case-sensitive; applied to the annotation name after the source prefix"))
    o.write("other\tany\t(no match, or matches of more than one category)\tambiguous matches are other\n")
CR = [(c, s, re.compile(rx)) for c, s, rx in CMAP]
def category(a):
    src, name = a.split(":", 1)
    hits = {c for c, s, rx in CR if (s == "any" and src != "COG" or s == src) and rx.search(name)}
    return hits.pop() if len(hits) == 1 else ("other (ambiguous)" if hits else "other")

# ---- per family output
bs = lambda x: bg_slot[x] / slots; bw3 = lambda x: bg_w3[x] / nmem
summ, conserved = {}, {}
oa = open(OUTS[1], "w")
oa.write("family_rep\tannotation\tn_animals\tanimals_pm3\tanimal_frac_pm3\tanimals_pm1\tanimal_frac_pm1\tn_members\tmembers_pm3\tmember_frac_pm3\tmembers_pm1\tmember_frac_pm1\tbg_slot_frac\tbg_win3_frac\tenrichment\tenrichment_win\tconserved50\tconserved75\tconserved90\tconserved50_win\tcategory\n")
for q in fams:
    nm, na, mnb, anb, a3m, a1m, a3a, a1a = fam_stats[q]; rows = []
    for x in sorted(a3a):
        af3 = a3a[x] / na; af1 = a1a.get(x, 0) / na; e = af3 / bs(x); ew = af3 / bw3(x)
        c = [af3 >= t and e >= 5 for t in (0.5, 0.75, 0.9)]; cw = af3 >= 0.5 and ew >= 5
        rows.append((x, af3, af1, e, ew, c, cw))
        oa.write("%s\t%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%d\t%.4f\t%d\t%.4f\t%.3g\t%.3g\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\n" % (
            q, x, na, a3a[x], af3, a1a.get(x, 0), af1, nm, a3m[x], a3m[x] / nm, a1m.get(x, 0), a1m.get(x, 0) / nm,
            bs(x), bw3(x), e, ew, *["yes" if v else "no" for v in c], "yes" if cw else "no", category(x)))
    cons = sorted([r for r in rows if r[5][0]], key=lambda r: -(r[1] * r[3]))
    conserved[q] = cons
    cats = [category(r[0]) for r in cons]
    mapped = [c for c in cats if not c.startswith("other")]
    prim = mapped[0] if mapped else ("other" if cons else "none")
    summ[q] = dict(nm=nm, na=na, mnb=mnb, anb=anb, c50=len(cons), c75=sum(r[5][1] for r in rows), c90=sum(r[5][2] for r in rows),
                   c50w=sum(r[6] for r in rows), prim=prim, cats=";".join(sorted(set(mapped))) or "-",
                   score=(cons[0][1] * cons[0][3]) if cons else 0.0,
                   top=";".join("%s(%.2f/%.2f/%.1f)" % (r[0], r[1], r[2], r[3]) for r in cons[:8]) or "-")
oa.close()

# ---- HIT families
best = {}
for db in ("PDB", "CATH50"):
    for l in open(R + "prostt5_main/eu15580_vs_%s.tsv" % db):
        f = l.rstrip("\n").split("\t"); e, b = float(f[2]), float(f[3])
        if f[0] not in best or (e, -b) < (best[f[0]][0], -best[f[0]][1]): best[f[0]] = (e, b, f[1], db)
hit = sorted(q for q, v in best.items() if v[0] <= 1e-3)
if len(hit) != 218: raise SystemExit("ABORT: expected 218 HIT families, got %d" % len(hit))
want = {best[q][2] for q in hit if best[q][3] == "PDB"}; title = {}
for e in open(FDB + "/PDB_h", "rb").read().decode("utf-8", "replace").split("\x00"):
    e = e.strip("\n")
    if not e: continue
    k, _, d = e.partition(" ")
    if k in want: title[k] = d.strip()
def desc(t, db):
    if db == "PDB": return title.get(t, "NA")
    m = re.search(r"_(\d+\.\d+\.\d+\.\d+)$", t)
    return ("CATH superfamily " + m.group(1)) if m else "CATH domain (no name in database)"
with open(OUTS[3], "w") as o:
    o.write("family_rep\ttop_target\tdatabase\tevalue\tbits\tdescription\tn_animals\tn_conserved50\tprimary_category\tconserved_annotations(animal_frac_pm3/animal_frac_pm1/enrichment)\n")
    for q in hit:
        e, b, t, db = best[q]; s = summ[q]
        o.write("%s\t%s\t%s\t%.3g\t%d\t%s\t%d\t%d\t%s\t%s\n" % (q, t, db, e, b, desc(t, db), s["na"], s["c50"], s["prim"],
                ";".join("%s(%.2f/%.2f/%.1f)" % (r[0], r[1], r[2], r[3]) for r in conserved[q]) or "-"))

with open(OUTS[0], "w") as o:
    o.write("family_rep\tn_members\tn_animals\tmembers_with_neighbour\tanimals_with_neighbour\tn_conserved50\tn_conserved75\tn_conserved90\tn_conserved50_win\tprimary_category\tcategories\tbest_score(animal_frac_pm3*enrichment)\tstructural_hit\ttop_conserved(animal_frac_pm3/animal_frac_pm1/enrichment)\n")
    for q in fams:
        s = summ[q]
        o.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\n" % (q, s["nm"], s["na"], s["mnb"], s["anb"], s["c50"], s["c75"],
                s["c90"], s["c50w"], s["prim"], s["cats"], s["score"], "yes" if q in set(hit) else "no", s["top"]))

# ---- report
print("\n== conserved neighbourhood (families with >= 1 conserved annotation)")
for k, lab in (("c50", ">=50% animals, enrichment >= 5 (slot background)"), ("c75", ">=75%"), ("c90", ">=90%"), ("c50w", ">=50%, enrichment_win >= 5 (window background)")):
    print("  %-52s %d" % (lab, sum(1 for s in summ.values() if s[k] > 0)))
print("\n== primary category (families with a conserved annotation at 50%)")
for c, n in collections.Counter(s["prim"] for s in summ.values()).most_common(): print("  %-32s %d" % (c, n))
print("== any category among conserved annotations")
cc = collections.Counter(c for s in summ.values() for c in s["cats"].split(";") if c != "-")
for c, n in cc.most_common(): print("  %-32s %d" % (c, n))
print("\n== top 20 families by animal_frac_pm3 x enrichment of their best conserved annotation")
for q in sorted(fams, key=lambda q: -summ[q]["score"])[:20]:
    s = summ[q]; print("  %-32s animals %3d  score %8.1f  %-28s %s" % (q, s["na"], s["score"], s["prim"], s["top"][:150]))
print("\n== HIT families: %d; with a conserved neighbourhood (50%%): %d" % (len(hit), sum(1 for q in hit if summ[q]["c50"])))
for q in hit:
    if summ[q]["c50"]:
        e, b, t, db = best[q]; print("  %-32s %-6s %-28s E %.2g  %-50s | %s" % (q, db, t, e, desc(t, db)[:50], summ[q]["top"][:120]))
print("\n== positive control %s" % POS)
s = summ[POS]; print("  members %d, animals %d, members with neighbour %d, animals with neighbour %d, conserved50 %d, primary %s" % (s["nm"], s["na"], s["mnb"], s["anb"], s["c50"], s["prim"]))
nm, na, mnb, anb, a3m, a1m, a3a, a1a = fam_stats[POS]
for x in ("Pfam:CRISPR_Cse1", "Pfam:Cas_CT1975", "Pfam:Cas_Cas5d"):
    print("  %-20s animals pm3 %d/%d pm1 %d/%d | members pm3 %d/%d pm1 %d/%d | enrichment %.1f, enrichment_win %.1f | category %s | conserved50 %s" % (
        x, a3a.get(x, 0), na, a1a.get(x, 0), na, a3m.get(x, 0), nm, a1m.get(x, 0), nm,
        (a3a.get(x, 0) / na) / bs(x) if bg_slot[x] else 0, (a3a.get(x, 0) / na) / bw3(x) if bg_w3[x] else 0, category(x),
        "yes" if any(r[0] == x for r in conserved[POS]) else "no"))
print("done")
