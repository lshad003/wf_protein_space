# Step 34f: final membership of the Cse2/CasB-like family and its type I-E cas loci (Figure 4d input).
# Usage: python3 step34f_cse2_loci.py <stage> <workdir>   stages run in order by scripts/step34f_cse2_loci.sh:
#   members  family members, per-contig strict_final check, contig extraction, list of contigs to map
#   verify   Prodigal (-p meta) proteins vs catalog proteins of the same name; aborts on any difference
#   loci     Pfam (project annotation via 95% reps + direct hmmsearch), MinCED arrays, per-contig table, ranking
# Inputs:
#   results/step34c2c_strict_v2_families.tsv   family UHM45.35758__k141_100597_236 (excluded == 0)
#   results/clusters_95_n88.tsv                95% rep (col 1) -> member (col 2)
#   results/step34a_rep_origin_v2_final.tsv    per supported rep: origin, fly95, frog95, fly90, lca_domain_v2,
#                                              strict_final. Origin is a property of the rep's own contig, so a member
#                                              contig takes the label of any supported rep on it (checked consistent).
#                                              strict_final at contig level = origin Tiara-prokaryotic AND fly90 == 0
#                                              AND lca_domain_v2 != Eukaryota (the strict_final rep column additionally
#                                              drops the 5 excluded families' reps, which are not on these contigs).
#   Contigs without a supported rep: the same rule from results/tiara/<stem>.tiara.tsv, a fresh minimap2 2.30
#     -x asm10 --secondary=no -c mapping to dm6 and the wood frog (coverage >= 0.50; identity >= 0.95 for fly95 and
#     frog95, >= 0.90 for fly90; as step34a), and the v2 lineage domain of the contig LCA (as step34c2c_taxv2.py).
#   metadata/samples_88_treatment.tsv (cohort), assemblies <FEC>/results/<stem>/scaffolds/<stem>_R.fa.gz,
#   catalog/db/LsPS_AA_all.fasta.gz (samtools faidx), results/pfam_hs/*.tblout (Pfam 38.2 hmmsearch --cut_ga of the
#   95% reps; source used by step29d), results/step29d_fold4_neighbours.tsv (the 37 members of Step 29d).
# Pfam: primary = the project annotation transferred from each gene's 95% rep (as step29d); check = hmmsearch --cut_ga
#   of the actual contig genes against the same Pfam-A 38.2 (hmmsearch, not hmmscan, per CLAUDE.md).
# Type I-E names: CRISPR_Cse1 cse1/casA; CRISPR_Cse2 cse2/casB; Cas_CT1975 cas7/casC; Cas_Cas5d cas5/casD;
#   CRISPR_assoc cas6e/casE; Cas_Cas1 cas1; CRISPR_Cas2 or Cas_Cas2CT1978 cas2; cas3 = Cas3_HD and/or DEAD plus
#   Helicase_C on one gene. The family member gene itself is the structural cse2 (Steps 25-26) and counts as cse2.
# Ranking: distinct type I-E genes present, then CRISPR array presence, then contig length. Genes are never combined
#   across contigs.
# Outputs: results/step34f_cse2_members.tsv, results/step34f_cse2_loci.tsv, results/fig_prep_cse2_locus_best.tsv,
#   and _rank2/_rank3 when no contig carries all 8 genes; intermediates in <workdir>.
import os, sys, re, gzip, glob, json, collections, subprocess
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
FEC = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal"
LCA = FEC + "/results_scaffold_classify_mmseqs"
REP = "UHM45.35758__k141_100597_236"
STAGE, D = sys.argv[1], sys.argv[2]
GRP = {"bacteria": "Tiara-prokaryotic", "archaea": "Tiara-prokaryotic", "prokarya": "Tiara-prokaryotic",
       "eukarya": "Tiara-eukaryotic", "organelle": "Tiara-eukaryotic", "unknown": "Tiara-unknown"}
def contig(g):
    s, rest = g.split("__", 1); return s, rest.rsplit("_", 1)[0]
def lineage_domain(lin):
    for t in lin.split(";"):
        n = t.split("_", 1)[-1]
        if n in ("Bacteria", "Archaea", "Eukaryota", "Viruses"): return n
    return "unclassified"
def paf_calls(p):
    iv, mt, bl, ql = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
    if os.path.exists(p):
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

if STAGE == "members":
    fam = [l.rstrip("\n").split("\t") for l in open(R + "step34c2c_strict_v2_families.tsv")]
    h = fam[0]; row = [dict(zip(h, r)) for r in fam[1:] if r[0] == REP]
    if len(row) != 1 or row[0]["excluded"] != "0": raise SystemExit("ABORT: family row")
    print("family row (step34c2c_strict_v2_families.tsv):", row[0])
    reps = {l.rstrip("\n").split("\t")[1] for l in open(R + "step34c_strict_eu_clusters30.tsv") if l.startswith(REP + "\t")}
    print("95% reps in the family (step34c_strict_eu_clusters30.tsv):", sorted(reps))
    mem = []
    with open(R + "clusters_95_n88.tsv") as fh:
        for l in fh:
            r, m = l.rstrip("\n").split("\t")
            if r in reps: mem.append(m)
    print("member proteins (clusters_95_n88.tsv):", len(mem), "| distinct:", len(set(mem)))
    mc = sorted({contig(m) for m in mem})
    # contig labels from supported reps on the same contig (step34a_rep_origin_v2_final.tsv)
    lab = collections.defaultdict(set)
    with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
        hh = next(fh).rstrip("\n").split("\t"); ix = {c: hh.index(c) for c in ("origin", "fly95", "frog95", "fly90", "lca_domain_v2")}
        want = set(mc)
        for l in fh:
            f = l.rstrip("\n").split("\t"); k = contig(f[0])
            if k in want: lab[k].add(tuple(f[ix[c]] for c in ("origin", "fly95", "frog95", "fly90", "lca_domain_v2")))
    incons = [k for k, v in lab.items() if len(v) > 1]
    if incons: raise SystemExit("ABORT: reps on one contig disagree: %s" % incons[:3])
    # Tiara and contig LCA for every member contig
    tia, tlen = {}, {}
    for s in sorted({k[0] for k in mc}):
        for i, l in enumerate(open(R + "tiara/%s.tiara.tsv" % s)):
            if i == 0: continue
            f = l.split("\t", 2); c = f[0].split()[0]
            if (s, c) in want: tia[(s, c)] = GRP[f[1]]; tlen[(s, c)] = int(re.search(r"len=(\d+)", f[0]).group(1))
    lca = {}
    for s in sorted({k[0] for k in mc}):
        p = "%s/%s/%s_uniref50_lca.tsv" % (LCA, s, s)
        if not os.path.exists(p): p += ".gz"
        if not os.path.exists(p): continue
        for l in (gzip.open(p, "rt") if p.endswith(".gz") else open(p)):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 10 and (s, f[0]) in want: lca[(s, f[0])] = (f[2], f[3], lineage_domain(f[9]))
    # extract contigs
    os.makedirs(D, exist_ok=True); seqlen = {}
    with open(D + "/contigs.fa", "w") as o:
        for s in sorted({k[0] for k in mc}):
            names = {c for ss, c in mc if ss == s}; w = False
            with gzip.open("%s/results/%s/scaffolds/%s_R.fa.gz" % (FEC, s, s), "rt") as fh:
                for l in fh:
                    if l[0] == ">":
                        c = l[1:].split()[0]; w = c in names
                        if w: o.write(">%s__%s\n" % (s, c)); seqlen[(s, c)] = 0
                    elif w: o.write(l); seqlen[(s, c)] += len(l.strip())
    miss = [k for k in mc if k not in seqlen]
    if miss: raise SystemExit("ABORT: contigs not found in assemblies: %s" % miss)
    tomap = [k for k in mc if k not in lab]
    with open(D + "/tomap.txt", "w") as o:
        for s, c in tomap: o.write("%s__%s\n" % (s, c))
    coh = {}
    for i, l in enumerate(open(W + "/metadata/samples_88_treatment.tsv")):
        if i: f = l.rstrip("\n").split("\t"); coh[f[0]] = f[1]
    json.dump({"mem": mem, "lab": {"%s__%s" % k: list(v)[0] for k, v in lab.items()},
               "tia": {"%s__%s" % k: v for k, v in tia.items()}, "lca": {"%s__%s" % k: v for k, v in lca.items()},
               "len": {"%s__%s" % k: v for k, v in seqlen.items()}, "coh": coh}, open(D + "/members.json", "w"))
    print("member contigs %d | with a supported-rep origin label %d | to map now %d | with Tiara call %d | with LCA row %d"
          % (len(mc), len(lab), len(tomap), len(tia), len(lca)))

elif STAGE == "verify":
    J = json.load(open(D + "/members.json"))
    pred = {}; h = None
    for l in open(D + "/prodigal.faa"):
        if l[0] == ">": h = l[1:].split()[0]; pred[h] = []
        else: pred[h].append(l.strip())
    pred = {k: "".join(v).rstrip("*") for k, v in pred.items()}
    with open(D + "/catalog_ids.txt", "w") as o:
        for k in sorted(pred): o.write(k + "\n")
    print("Prodigal proteins on member contigs:", len(pred))

elif STAGE == "compare":
    pred = {}; h = None
    for l in open(D + "/prodigal.faa"):
        if l[0] == ">": h = l[1:].split()[0]; pred[h] = []
        else: pred[h].append(l.strip())
    pred = {k: "".join(v).rstrip("*") for k, v in pred.items()}
    cat = {}; h = None
    for l in open(D + "/catalog.faa"):
        if l[0] == ">": h = l[1:].split()[0]; cat[h] = []
        else: cat[h].append(l.strip())
    cat = {k: "".join(v).rstrip("*") for k, v in cat.items()}
    missing = [k for k in pred if k not in cat]; diff = [k for k in pred if k in cat and cat[k] != pred[k]]
    print("compared %d | missing from catalog %d | sequence differs %d" % (len(pred), len(missing), len(diff)))
    if missing or diff:
        print("missing:", missing[:10]); print("differ:", diff[:10]); raise SystemExit("ABORT: Prodigal proteins differ from the catalog")
    J = json.load(open(D + "/members.json"))
    bad = [m for m in J["mem"] if m not in pred]
    if bad: raise SystemExit("ABORT: family members not predicted: %s" % bad)
    print("all %d family members predicted and identical" % len(J["mem"]))

elif STAGE == "loci":
    J = json.load(open(D + "/members.json"))
    mem = J["mem"]; memset = set(mem)
    # genes from Prodigal GFF
    genes = collections.defaultdict(list)
    for l in open(D + "/prodigal.gff"):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        n = re.search(r"ID=\d+_(\d+)", f[8]).group(1); pt = re.search(r"partial=(\d\d)", f[8]).group(1)
        genes[f[0]].append(dict(gene="%s_%s" % (f[0], n), start=int(f[3]), end=int(f[4]), strand=f[6], partial=pt))
    allg = {g["gene"] for v in genes.values() for g in v}
    # project Pfam via 95% reps
    rep95 = {}
    with open(R + "clusters_95_n88.tsv") as fh:
        for l in fh:
            r, m = l.rstrip("\n").split("\t")
            if m in allg: rep95[m] = r
    reps = set(rep95.values()); pf = collections.defaultdict(set)
    for f in glob.glob(R + "pfam_hs/*.tblout"):
        for l in open(f):
            if l[0] == "#": continue
            p = l.split()
            if p[0] in reps: pf[p[0]].add(p[2])
    direct = collections.defaultdict(set)
    for l in open(D + "/hmmsearch.tblout"):
        if l[0] == "#": continue
        p = l.split(); direct[p[0]].add(p[2])
    print("genes on member contigs %d | with a 95%% rep %d | rep with Pfam %d | direct hmmsearch with Pfam %d"
          % (len(allg), len(rep95), sum(1 for g in allg if pf.get(rep95.get(g))), sum(1 for g in allg if direct.get(g))))
    agree = sum(1 for g in allg if pf.get(rep95.get(g), set()) == direct.get(g, set()))
    print("genes where rep-transferred and direct Pfam sets agree: %d of %d" % (agree, len(allg)))
    names_seen = collections.Counter(x for g in allg for x in (pf.get(rep95.get(g), set()) | direct.get(g, set()))
                                     if re.search(r"cas|crispr|cse|DEAD|Helicase_C|HD", x, re.I))
    print("cas/CRISPR/helicase/HD Pfam names found (rep or direct):", dict(names_seen))
    MAP = [("cse1", {"CRISPR_Cse1"}), ("cse2", {"CRISPR_Cse2"}), ("cas7", {"Cas_CT1975"}), ("cas5", {"Cas_Cas5d"}),
           ("cas6e", {"CRISPR_assoc"}), ("cas1", {"Cas_Cas1"}), ("cas2", {"CRISPR_Cas2", "Cas_Cas2CT1978"})]
    def casname(g, d):
        n = [c for c, s in MAP if d & s]
        if "Cas3_HD" in d or ("DEAD" in d and "Helicase_C" in d): n.append("cas3")
        if g in memset and "cse2" not in n: n.append("cse2")
        return n
    # MinCED arrays
    arr = collections.defaultdict(list)
    if os.path.exists(D + "/minced.gff"):
        for l in open(D + "/minced.gff"):
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] == "repeat_region":
                m = re.search(r"score=?(\d+)", f[8]); sc = f[5]
                arr[f[0]].append(dict(start=int(f[3]), end=int(f[4]), n_repeats=int(sc) if sc.isdigit() else None))
    # repeat counts from the text report (authoritative)
    if os.path.exists(D + "/minced.txt"):
        cur = None; cnt = collections.defaultdict(list)
        for l in open(D + "/minced.txt"):
            m = re.match(r"Sequence '(\S+)'", l)
            if m: cur = m.group(1)
            m = re.match(r"CRISPR (\d+)\s+Range: (\d+) - (\d+)", l)
            if m: cnt[cur].append([int(m.group(2)), int(m.group(3)), None])
            m = re.match(r"Repeats: (\d+)", l.strip())
            if m and cnt[cur]: cnt[cur][-1][2] = int(m.group(1))
        for c, v in cnt.items():
            arr[c] = [dict(start=a, end=b, n_repeats=n) for a, b, n in v]
    # Origin per contig and table
    tomap = {}
    for p, k in ((D + "/fly.paf", "fly"), (D + "/frog.paf", "frog")):
        for c, (cv, idn) in paf_calls(p).items(): tomap.setdefault(c, {})[k] = (cv, idn)
    rows = []; IE = ["cse1", "cse2", "cas7", "cas5", "cas6e", "cas1", "cas2", "cas3"]
    stems = collections.Counter(m.split("__")[0] for m in mem)
    for m in sorted(mem):
        s, c = contig(m); key = "%s__%s" % (s, c)
        if key in J["lab"]:
            o, f95, fr95, f90, dom = J["lab"][key]; src = "step34a_rep_origin_v2_final"
            strict = o == "Tiara-prokaryotic" and f90 == "0" and dom != "Eukaryota"
        else:
            t = J["tia"].get(key, "<1 kb unclassified"); mp = tomap.get(key, {})
            fl = mp.get("fly", (0, 0)); fr = mp.get("frog", (0, 0))
            f95 = str(int(fl[0] >= 0.5 and fl[1] >= 0.95)); f90 = str(int(fl[0] >= 0.5 and fl[1] >= 0.90)); fr95 = str(int(fr[0] >= 0.5 and fr[1] >= 0.95))
            dom = J["lca"].get(key, [None, None, "NA"])[2]
            o = "fly" if f95 == "1" else ("frog" if fr95 == "1" else t)
            strict = o == "Tiara-prokaryotic" and f90 == "0" and dom != "Eukaryota"; src = "computed (Tiara + minimap2 -c + LCA v2)"
        gl = sorted(genes.get(key, []), key=lambda g: g["start"])
        ann = []
        for i, g in enumerate(gl, 1):
            d = pf.get(rep95.get(g["gene"]), set()); dd = direct.get(g["gene"], set())
            ann.append(dict(g, idx=i, pfam=sorted(d), pfam_direct=sorted(dd), cas=casname(g["gene"], d | dd)))
        pos = [a["idx"] for a in ann if a["gene"] == m][0]
        present = sorted({n for a in ann for n in a["cas"]}, key=IE.index)
        casg = [a for a in ann if a["cas"]]
        edge = [a["gene"].rsplit("_", 1)[1] for a in casg if a["partial"] != "00"]
        A = arr.get(key, [])
        rows.append(dict(member=m, stem=s, cohort=J["coh"].get(s, "NA"), animal=s.split(".")[0], contig=c, contig_len=J["len"][key],
                         n_genes=len(gl), cse2_gene_index=pos, n_typeIE=len(present), typeIE_present=",".join(present) or "-",
                         typeIE_missing=",".join(x for x in IE if x not in present) or "-",
                         cas_genes_partial=",".join(edge) or "-", first_gene_partial=gl[0]["partial"] if gl else "NA",
                         last_gene_partial=gl[-1]["partial"] if gl else "NA", crispr_array=int(bool(A)),
                         arrays=";".join("%d-%d:%s" % (a["start"], a["end"], a["n_repeats"]) for a in A) or "-",
                         origin=o, fly95=f95, fly90=f90, frog95=fr95, lca_domain_v2=dom, origin_source=src,
                         strict_final_contig=int(strict), _ann=ann, _arr=A))
    keep = [r for r in rows if r["strict_final_contig"]]
    print("\nmembers %d | on strict_final contigs %d" % (len(rows), len(keep)))
    for r in rows:
        if not r["strict_final_contig"]: print("  not strict_final:", r["member"], r["origin"], r["fly90"], r["lca_domain_v2"], r["origin_source"])
    keep.sort(key=lambda r: (-r["n_typeIE"], -r["crispr_array"], -r["contig_len"], r["member"]))
    for i, r in enumerate(keep, 1): r["rank"] = i
    C = ["rank", "member", "stem", "cohort", "animal", "contig", "contig_len", "n_genes", "cse2_gene_index", "n_typeIE",
         "typeIE_present", "typeIE_missing", "cas_genes_partial", "first_gene_partial", "last_gene_partial", "crispr_array",
         "arrays", "origin", "fly95", "fly90", "frog95", "lca_domain_v2", "origin_source", "strict_final_contig"]
    with open(R + "step34f_cse2_loci.tsv", "w") as o:
        o.write("\t".join(C) + "\n")
        for r in keep: o.write("\t".join(str(r[c]) for c in C) + "\n")
    with open(R + "step34f_cse2_members.tsv", "w") as o:
        o.write("member\tstem\tcohort\tanimal\tcontig\tcontig_len\torigin\tfly95\tfly90\tfrog95\tlca_domain_v2\torigin_source\tstrict_final_contig\n")
        for r in rows: o.write("\t".join(str(r[c]) for c in ("member", "stem", "cohort", "animal", "contig", "contig_len", "origin", "fly95", "fly90", "frog95", "lca_domain_v2", "origin_source", "strict_final_contig")) + "\n")
    print("strict_final members %d | metagenomes %d | animals %d | cohorts %s"
          % (len(keep), len({r["stem"] for r in keep}), len({r["animal"] for r in keep}), dict(collections.Counter(r["cohort"] for r in keep))))
    s29 = {l.split("\t")[2] for i, l in enumerate(open(R + "step29d_fold4_neighbours.tsv")) if i and l.startswith("UHM636.41030__k141_332084_5\t")}
    print("Step 29d members: %d | in the 45: %d | of the 45 not in Step 29d: %d" % (len(s29), len(s29 & memset), len(memset - s29)))
    for m in sorted(memset - s29):
        a = [r for r in rows if r["member"] == m][0]; g = [x for x in a["_ann"] if x["gene"] == m][0]
        print("   not in 29d:", m, "partial", g["partial"], "length_aa", (g["end"] - g["start"] + 1) // 3 - (1 if g["partial"][1] == "0" else 0))
    print("\ntop contigs:")
    for r in keep[:5]:
        print("  #%d %s__%s len %d genes %d cse2 at %d | %d type I-E: %s | missing %s | partial cas genes %s | array %s"
              % (r["rank"], r["stem"], r["contig"], r["contig_len"], r["n_genes"], r["cse2_gene_index"], r["n_typeIE"],
                 r["typeIE_present"], r["typeIE_missing"], r["cas_genes_partial"], r["arrays"]))
    def write(r, p):
        with open(p, "w") as o:
            o.write("type\tgene\tstart\tend\tstrand\tpfam\tpfam_direct\tcas_name\tpartial\tn_repeats\n")
            for a in r["_ann"]:
                o.write("gene\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t-\n" % (a["gene"], a["start"], a["end"], a["strand"], ",".join(a["pfam"]) or "-",
                        ",".join(a["pfam_direct"]) or "-", ",".join(a["cas"]) or "-", a["partial"]))
            for i, x in enumerate(r["_arr"], 1):
                o.write("array\tCRISPR%d\t%d\t%d\t.\t-\t-\t-\t-\t%s\n" % (i, x["start"], x["end"], x["n_repeats"]))
        print("wrote", p)
    write(keep[0], R + "fig_prep_cse2_locus_best.tsv")
    if keep[0]["n_typeIE"] < 8:
        for k in (2, 3):
            if len(keep) >= k: write(keep[k - 1], R + "fig_prep_cse2_locus_rank%d.tsv" % k)
    print("wrote", R + "step34f_cse2_loci.tsv", "and", R + "step34f_cse2_members.tsv")
