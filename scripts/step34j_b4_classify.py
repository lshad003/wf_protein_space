# Step 34j B4: structure hit classes for the ESMFold models of the 4,250 member-validated EU families and the 4,509
#   known controls (results/step34j/esmfold_input_labels.tsv).
# Model confidence: mean C-alpha B-factor of the model PDB (ESMFold pLDDT, 0-100); pTM and ESMFold's own mean pLDDT from
#   the ESMFold logs (logs/step34j_b1_esmfold_cpu.*.log, "Predicted structure for ... pLDDT x, pTM y").
# Hits: results/step34j/foldseek/{PDB,CATH50,AFDB50,ESMAtlas30}.m8 (scripts/step34j_b3_foldseek.sh), full models and
#   Merizo domains (<model>__D<k>); a hit is E <= 1e-3. A model's hits are the union of its own and its domains' hits.
# AFDB annotation: the 25 best AFDB50 hits (by E) per query unit; UniProt entries fetched in batches of 500 with curl from
#   the UniProt REST API (accessions endpoint; fields protein_name, xref_pfam, xref_interpro; raw responses kept in
#   results/step34j/uniprot/). Annotated = informative protein name, or >= 1 Pfam not named DUF/UPF/"unknown function"
#   (Pfam 38.2 Pfam-A.hmm.dat), or >= 1 InterPro entry not named unknown function/uncharacterised (InterPro entry.list,
#   results/step34j/ref_interpro). Uninformative name: uncharacteri[sz]ed, hypothetical, unknown function, DUF<n>,
#   UPF<n>, "predicted protein", "putative protein", or a name made only of a DUF/UPF domain-containing phrase.
# Hit class (precedence): experimental (PDB or CATH50) > AFDB annotated > AFDB unannotated only > ESM Atlas only > none.
# Structure-confident annotated hit: class experimental or AFDB annotated, with a hit of that class at qTM >= 0.5 and
#   prob >= 0.5 (Step 29e pass rule) and model mean C-alpha pLDDT >= 70.
# Names: PDB entry COMPOUND from wwPDB entries.idx; CATH superfamily name from cath-names.txt (v4.4.0; CATH domain IDs
#   via cath-domain-list.txt, af_ targets carry the superfamily code).
# ProstT5 comparison: EU = results/step34d_prostt5/queries_vs_{PDB,CATH50}.tsv; controls = Step 34d known rows of the
#   fixed set plus results/step34d2_prostt5/krepl_vs_*.tsv; hit = best E <= 1e-3.
# Outputs (refuse to overwrite): results/step34j_structures.tsv (one row per model), results/step34j/structures_by_bin.tsv,
#   results/step34j/structures_vs_prostt5.tsv, results/step34j/uniprot/.
# 2026-10-22 patch (UniProt gap; run 1 = job 30171185 wrote the outputs above with 31,660 of 88,365 accessions resolved):
#   no download; reads the saved responses in results/step34j/uniprot/ (aborts unless the batch count matches the
#   requested list), maps secondary accessions back to the requested ones (results/step34j/uniprot_unisave/sec_acc/*.json.gz,
#   scripts/step34j_b4a_unisave_fetch.py part A), and for accessions still unresolved reads the last UniSave version
#   (uniprot_unisave/unisave_raw.jsonl, part B; name = first top-level DE RecName/SubName Full=, Pfam/InterPro from DR
#   lines). Same annotation rule for all three sources. New column afdb_annot_sources (sources of the annotated top-25
#   accessions: live, secondary, unisave). Outputs under new names (_v2); run 1 outputs are kept.
# 2026-10-23 ProtNLM rule: a name whose evidence is only ECO:0008006 (Google ProtNLM, predicted) does not count as an
#   informative name; the entry is then annotated only through Pfam/InterPro. Live entries: evidence of the first
#   recommendedName/submissionNames fullName in the saved full JSON (sec_acc/*.json.gz, which holds every live entry);
#   secondary: that of the primary; UniSave: the {...} evidence of the DE name line. hit_class, structure_confident_annotated
#   and the confident/best columns follow this rule (c1 uses them); hit_class_protnlm_annotated and
#   structure_confident_protnlm_annotated give the rule that accepts ProtNLM names. structures_by_bin_v2.tsv has a column
#   protnlm_rule (unannotated / annotated) with both.
import os, re, glob, subprocess, collections, time, json
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34j/"
OS, OB, OP, UD = R + "step34j_structures_v2.tsv", D + "structures_by_bin_v2.tsv", D + "structures_vs_prostt5_v2.tsv", D + "uniprot/"
US = D + "uniprot_unisave/"
for p in (OS, OB, OP, D + "uniprot_annotation_v2.tsv"):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
for p in (UD, US + "sec_acc", US + "unisave_raw.jsonl"):
    if not os.path.exists(p): raise SystemExit("ABORT: missing input " + p)
E_HIT = 1e-3
lab = {}
with open(D + "esmfold_input_labels.tsv") as fh:
    next(fh)
    for l in fh:
        i, s, fam, L, _ = l.rstrip("\n").split("\t"); lab[i] = (s, fam, int(L))
print("models expected %d" % len(lab))
plddt, nres = {}, {}
for i in lab:
    p = D + "esmfold_pdb/%s.pdb" % i
    if not os.path.exists(p): continue
    b = [float(l[60:66]) for l in open(p) if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    plddt[i] = sum(b) / len(b); nres[i] = len(b)
mx = max(plddt.values()); print("models found %d; max mean CA B-factor %.2f" % (len(plddt), mx))
if mx <= 1.0: plddt = {k: v * 100 for k, v in plddt.items()}   # guard: 0-1 scale
bad = [i for i in plddt if nres[i] != lab[i][2]]
if bad: raise SystemExit("ABORT: %d models with residue count != sequence length" % len(bad))
ptm, eplddt = {}, {}
rx = re.compile(r"Predicted structure for (\S+) with length (\d+), pLDDT ([\d.]+), pTM ([\d.]+)")
for f in glob.glob(W + "/logs/step34j_b1_esmfold_cpu.*.log") + glob.glob(W + "/logs/step34j_b1_esmfold_rerun*.log"):
    for l in open(f):
        m = rx.search(l)
        if m and m.group(1) in lab: ptm[m.group(1)] = float(m.group(4)); eplddt[m.group(1)] = float(m.group(3))
print("pTM parsed for %d models" % len(ptm))
units = set(plddt) | {os.path.basename(x)[:-4] for x in glob.glob(D + "domains/*.pdb")}
def unit(q):
    for c in (q, re.sub(r"\.pdb$", "", q), re.sub(r"(\.pdb)?_[A-Z]$", "", q)):
        if c in units: return c
    raise SystemExit("ABORT: cannot map query %s" % q)
def model_of(u): return u.split("__D")[0] if re.search(r"__D\d+$", u) else u
hits = {db: collections.defaultdict(list) for db in ("PDB", "CATH50", "AFDB50", "ESMAtlas30")}
for db in hits:
    p = D + "foldseek/%s.m8" % db
    if not os.path.exists(p): raise SystemExit("ABORT: missing " + p)
    n = 0
    for l in open(p):
        f = l.rstrip("\n").split("\t")
        e = float(f[6])
        if e > E_HIT: continue
        u = unit(f[0]); n += 1
        hits[db][model_of(u)].append(dict(unit=u, target=f[1], e=e, bits=float(f[7]), prob=float(f[8]), qtm=float(f[10]), ttm=float(f[11]), qcov=float(f[13])))
    print("%s: hit rows (E <= 1e-3) %d, models with a hit %d" % (db, n, len(hits[db])))
# AFDB accessions
def acc(t):
    m = re.match(r"AF-([A-Z0-9]+)-F\d+", t); return m.group(1) if m else None
want = set()
for m, hs in hits["AFDB50"].items():
    byu = collections.defaultdict(list)
    for h in hs: byu[h["unit"]].append(h)
    for u, v in byu.items():
        for h in sorted(v, key=lambda h: (h["e"], -h["bits"]))[:25]:
            a = acc(h["target"]); h["top25"] = True
            if a: want.add(a)
nacc_fail = sum(1 for hs in hits["AFDB50"].values() for h in hs if h.get("top25") and not acc(h["target"]))
print("AFDB accessions to fetch %d (top-25 hits whose target name gave no accession: %d)" % (len(want), nacc_fail))
want = sorted(want); up = {}; src = {}
nb = (len(want) + 499) // 500
if len(glob.glob(UD + "batch_*.tsv")) != nb: raise SystemExit("ABORT: %d saved batches, %d expected" % (len(glob.glob(UD + "batch_*.tsv")), nb))
live = {}
for bi in range(0, len(want), 500):
    out = UD + "batch_%05d.tsv" % (bi // 500)
    for l in list(open(out))[1:]:
        f = l.rstrip("\n").split("\t") + ["", "", ""]
        live[f[0]] = (f[1], [x for x in f[2].split(";") if x], [x for x in f[3].split(";") if x])
for a in want:
    if a in live: up[a] = live[a]; src[a] = "live"
print("UniProt saved entries %d; requested accessions answered directly %d of %d" % (len(live), len(up), len(want)))
sec = {}; pnlm_live = {}; pnlm = {}
import gzip
for f in sorted(glob.glob(US + "sec_acc/batch_*.json.gz")):
    for x in json.load(gzip.open(f, "rt"))["results"]:
        for s2 in x.get("secondaryAccessions", []): sec.setdefault(s2, set()).add(x["primaryAccession"])
        pd = x.get("proteinDescription", {})
        fn = pd["recommendedName"]["fullName"] if "recommendedName" in pd else pd["submissionNames"][0]["fullName"] if pd.get("submissionNames") else {}
        ev = [e.get("evidenceCode") for e in fn.get("evidences", [])]
        pnlm_live[x["primaryAccession"]] = bool(ev) and all(e == "ECO:0008006" for e in ev)
nsec = 0
for a in want:
    if a in up or a not in sec: continue
    p = [q for q in sec[a] if q in live]
    if len(p) == 1: up[a] = live[p[0]]; src[a] = "secondary"; nsec += 1; pnlm[a] = pnlm_live[p[0]]
for a in want:
    if src.get(a) == "live": pnlm[a] = pnlm_live[a]
print("mapped back via secondary accession %d (saved entries with sec_acc lists %d)" % (nsec, len({q for v in sec.values() for q in v})))
def flat(t):
    nm = None; ev = ""; pf, ip = [], []; sub = False
    for l in t.splitlines():
        if l.startswith("DE   "):
            s2 = l[5:]
            if s2.startswith(("Contains:", "Includes:")): sub = True
            m = re.match(r"\s*(RecName|SubName): Full=(.*?)(\s*\{[^}]*\})?;\s*$", s2)
            if m and not sub and nm is None and not s2.startswith(" "): nm = m.group(2); ev = m.group(3) or ""
        elif l.startswith("DR   Pfam; "): pf.append(l.split(";")[1].strip())
        elif l.startswith("DR   InterPro; "): ip.append(l.split(";")[1].strip())
    codes = re.findall(r"ECO:\d+", ev)
    return nm, pf, ip, bool(codes) and all(c == "ECO:0008006" for c in codes)
uns = collections.Counter(); unsrec = {}
for l in open(US + "unisave_raw.jsonl"):
    try: r = json.loads(l)
    except ValueError: continue   # a line cut by a cancelled fetch run; that accession is refetched on resume
    if r["status"] != "failed" or r["acc"] not in unsrec: unsrec[r["acc"]] = r
nus = 0
for a in want:
    if a in up: continue
    r = unsrec.get(a)
    if r is None: uns["absent"] += 1; continue
    uns[r["status"]] += 1
    if r["status"] != "ok": continue
    nm, pf, ip, pn = flat(r["txt_raw"])
    up[a] = (nm or "", pf, ip); src[a] = "unisave"; nus += 1; pnlm[a] = pn
print("UniSave for the still-unresolved: %s; resolved %d" % (dict(uns), nus))
print("UniProt entries resolved %d of %d requested (live %d, secondary %d, unisave %d)" % (len(up), len(want),
      sum(1 for v in src.values() if v == "live"), nsec, nus))
pfn = {}; cur = None
for l in open("/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat"):
    if l.startswith("#=GF ID"): cur = [l.split()[2], None]
    elif l.startswith("#=GF AC"): pfn[l.split()[2].split(".")[0]] = cur
    elif l.startswith("#=GF DE"): cur[1] = l[10:].strip()
ipn = {}
for l in list(open(D + "ref_interpro/entry.list"))[1:]:
    f = l.rstrip("\n").split("\t"); ipn[f[0]] = f[2]
UNK = re.compile(r"uncharacteri[sz]ed|hypothetical|unknown function|\bDUF\d|\bUPF\d|predicted protein|putative protein|^protein$|^$", re.I)
def informative_name(n):
    n = re.sub(r"\s*\(.*", "", n).strip()
    return not UNK.search(n)
def pf_inf(a):
    x = pfn.get(a)
    return bool(x) and not re.match(r"(DUF|UPF)\d", x[0]) and not re.search(r"unknown function", x[1] or "", re.I)
def ip_inf(a):
    x = ipn.get(a)
    return bool(x) and not re.search(r"unknown function|uncharacteri[sz]ed|\bDUF\d|\bUPF\d", x, re.I)
def annotated(a, strict=True):
    if a not in up: return False
    n, pf, ip = up[a]
    return (informative_name(n) and not (strict and pnlm[a])) or any(pf_inf(x) for x in pf) or any(ip_inf(x) for x in ip)
print("\nresolved entries by source: n | informative name | of those ProtNLM-only | non-ProtNLM informative name | informative Pfam/InterPro | annotated (ProtNLM as unannotated) | annotated (ProtNLM accepted)")
for sname in ("live", "secondary", "unisave", "all"):
    A = [a for a in up if sname == "all" or src[a] == sname]
    inf = [a for a in A if informative_name(up[a][0])]
    print("  %-9s %6d | %6d | %6d | %6d | %6d | %6d | %6d" % (sname, len(A), len(inf), sum(pnlm[a] for a in inf), sum(not pnlm[a] for a in inf),
          sum(any(pf_inf(x) for x in up[a][1]) or any(ip_inf(x) for x in up[a][2]) for a in A), sum(annotated(a) for a in A), sum(annotated(a, False) for a in A)))
with open(D + "uniprot_annotation_v2.tsv", "w") as o:
    o.write("accession\tsource\tname\tname_informative\tname_protnlm_only\tpfam\tinterpro\tannotated_protnlm_unannotated\tannotated_protnlm_accepted\n")
    for a in sorted(up):
        n, pf, ip = up[a]
        o.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\n" % (a, src[a], n.replace("\t", " ") or "-", informative_name(n), pnlm[a], ";".join(pf) or "-", ";".join(ip) or "-", annotated(a), annotated(a, False)))
pdbt = {}
for l in open(D + "ref_pdb/entries.idx", encoding="latin-1"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 4 and len(f[0]) == 4: pdbt[f[0].upper()] = f[3]
cdom, cname = {}, {}
for l in open(D + "ref_cath/cath-domain-list.txt"):
    if l[0] == "#": continue
    f = l.split(); cdom[f[0]] = ".".join(f[1:5])
for l in open(D + "ref_cath/cath-names.txt"):
    if l[0] == "#": continue
    f = l.rstrip("\n").split(None, 2)
    if len(f) == 3: cname[f[0]] = f[2].lstrip(":")
def tname(db, t):
    if db == "PDB": return pdbt.get(t[:4].upper(), "NA")
    if db == "CATH50":
        code = t.rsplit("_", 1)[1] if t.startswith("af_") else cdom.get(t, "NA")
        return "%s %s" % (code, cname.get(code, "NA"))
    if db == "AFDB50":
        a = acc(t); return up[a][0] if a in up else "NA"
    return "-"
def best(v): return min(v, key=lambda h: (h["e"], -h["bits"])) if v else None
def ok(h): return h["qtm"] >= 0.5 and h["prob"] >= 0.5
prost = set()
for f in ("queries_vs_PDB.tsv", "queries_vs_CATH50.tsv", "known_vs_PDB.tsv", "known_vs_CATH50.tsv"):
    for l in open(R + "step34d_prostt5/" + f):
        x = l.split("\t")
        if float(x[2]) <= E_HIT: prost.add(x[0])
for f in ("krepl_vs_PDB.tsv", "krepl_vs_CATH50.tsv"):
    for l in open(R + "step34d2_prostt5/" + f):
        x = l.split("\t")
        if float(x[2]) <= E_HIT: prost.add(x[0])
def lbin(L): return "<100" if L < 100 else "100-199" if L < 200 else "200-399" if L < 400 else ">=400"
def pbin(p): return "NA" if p is None else "<50" if p < 50 else "50-70" if p < 70 else ">=70"
ndom = collections.Counter(model_of(u) for u in units if "__D" in u)
rows = []
for i, (s, fam, L) in sorted(lab.items()):
    p = plddt.get(i)
    ex = hits["PDB"].get(i, []) + hits["CATH50"].get(i, [])
    af = hits["AFDB50"].get(i, []); es = hits["ESMAtlas30"].get(i, [])
    def classify(strict):
        afa = [h for h in af if h.get("top25") and annotated(acc(h["target"]) or "", strict)]
        if p is None: cls = "no_model"
        elif ex: cls = "experimental"
        elif afa: cls = "AFDB_annotated"
        elif af: cls = "AFDB_unannotated_only"
        elif es: cls = "ESMAtlas_only"
        else: cls = "none"
        top = {"experimental": ex, "AFDB_annotated": afa, "AFDB_unannotated_only": af, "ESMAtlas_only": es}.get(cls, [])
        conf = int(cls in ("experimental", "AFDB_annotated") and p is not None and p >= 70 and any(ok(h) for h in top))
        return afa, cls, top, conf
    afa, cls, top, conf = classify(True)
    _, cls2, _, conf2 = classify(False)
    b = best(top)
    def dbcols(db):
        h = best(hits[db].get(i, []))
        return [h["target"], "%.2e" % h["e"], "%.3f" % h["qtm"], "%.3f" % h["prob"], h["unit"]] if h else ["-"] * 5
    bdb = None
    if b: bdb = "PDB" if b in hits["PDB"].get(i, []) else "CATH50" if b in hits["CATH50"].get(i, []) else "AFDB50" if b in af else "ESMAtlas30"
    bc = max((h for h in top if ok(h)), key=lambda h: h["bits"], default=None) if conf else None
    bcdb = None
    if bc: bcdb = "PDB" if bc in hits["PDB"].get(i, []) else "CATH50" if bc in hits["CATH50"].get(i, []) else "AFDB50"
    rows.append([i, s, fam, L, lbin(L), "%.2f" % p if p is not None else "NA", eplddt.get(i, "NA"), ptm.get(i, "NA"), pbin(p), ndom.get(i, 0),
                 cls, bdb or "-", b["target"] if b else "-", "%.2e" % b["e"] if b else "-", tname(bdb, b["target"]) if b else "-",
                 conf, bcdb or "-", bc["target"] if bc else "-", "%.3f" % bc["qtm"] if bc else "-", tname(bcdb, bc["target"]) if bc else "-",
                 int(i in prost)] + dbcols("PDB") + dbcols("CATH50") + dbcols("AFDB50") + dbcols("ESMAtlas30")
                + [",".join(sorted({src[acc(h["target"])] for h in afa})) or "-", cls2, conf2])
H = ["id", "set", "family_rep", "length_aa", "length_bin", "mean_ca_plddt", "esmfold_plddt", "ptm", "plddt_bin", "n_merizo_domains",
     "hit_class", "best_db", "best_target", "best_evalue", "best_target_name", "structure_confident_annotated", "confident_db",
     "confident_target", "confident_qtm", "confident_target_name", "prostt5_hit"]
for db in ("PDB", "CATH50", "AFDB50", "ESMAtlas30"): H += ["%s_%s" % (db, x) for x in ("target", "evalue", "qtm", "prob", "unit")]
H += ["afdb_annot_sources", "hit_class_protnlm_annotated", "structure_confident_protnlm_annotated"]
with open(OS, "w") as o:
    o.write("\t".join(H) + "\n")
    for r in rows: o.write("\t".join(map(str, r)).replace("\n", " ") + "\n")
CL = ["experimental", "AFDB_annotated", "AFDB_unannotated_only", "ESMAtlas_only", "none", "no_model"]
with open(OB, "w") as o:
    o.write("protnlm_rule\tset\tbin_type\tbin\tn\t" + "\t".join("%s\tpct_%s" % (c, c) for c in CL) + "\tstructure_confident_annotated\n")
    for rule, ci, fi in (("unannotated", 10, 15), ("annotated", -2, -1)):
        for s in ("EU", "control"):
            for bt, idx, bins in (("all", None, ["all"]), ("length", 4, ["<100", "100-199", "200-399", ">=400"]), ("plddt", 8, ["<50", "50-70", ">=70"])):
                for bn in bins:
                    R2 = [r for r in rows if r[1] == s and (idx is None or r[idx] == bn)]
                    c = collections.Counter(r[ci] for r in R2); n = len(R2)
                    o.write("%s\t%s\t%s\t%s\t%d\t%s\t%d\n" % (rule, s, bt, bn, n, "\t".join("%d\t%.2f" % (c[k], 100.0 * c[k] / n if n else 0) for k in CL), sum(r[fi] for r in R2)))
                    print("ProtNLM %-11s %-7s %-6s %-8s n=%-5d %s | confident %d" % (rule, s, bt, bn, n, " ".join("%s %d" % (k, c[k]) for k in CL), sum(r[fi] for r in R2)))
with open(OP, "w") as o:
    o.write("set\tprostt5_hit\tn\t" + "\t".join(CL) + "\n")
    for s in ("EU", "control"):
        for ph in (1, 0):
            R2 = [r for r in rows if r[1] == s and r[20] == ph]; c = collections.Counter(r[10] for r in R2)
            o.write("%s\t%d\t%d\t%s\n" % (s, ph, len(R2), "\t".join(str(c[k]) for k in CL)))
            print("ProstT5 %s hit=%d n=%d: %s" % (s, ph, len(R2), " ".join("%s %d" % (k, c[k]) for k in CL)))
print("wrote", OS, OB, OP, D + "uniprot_annotation_v2.tsv")
