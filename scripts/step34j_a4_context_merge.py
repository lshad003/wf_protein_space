# Step 34j A4: Part A table. Per recurrent EU family: NR contexts as specified (95%-rep lists, >= 80% shared;
#   results/step34j/context_scores.tsv, column n_nr_contexts_rep95, reported only), and the scores used for calls over
#   containment-collapsed backgrounds on 50%-cluster lists (results/step34j/context_scores_nr50c.tsv, A2c; the 95%-rep
#   collapse saturated, see step34j_a2b/a2c headers): top-neighbour fraction with its control percentile (saturated,
#   reported), mean pairwise Jaccard of neighbour 50%-family sets with its control percentile (context_jaccard.tsv, A2d).
#   Conserved context (criterion): >= 3 backgrounds AND Jaccard percentile vs controls >= 90. Theme
#   fractions over backgrounds (seeds: results/step34j/context_nr50c_seeds.tsv; neighbour annotation: context_neighbours.tsv,
#   rep_annot.tsv; themes: scripts/step34j_themes.py), geNomad and DefenseFinder status of member contigs, and for the
#   186 ProstT5/Foldseek structural-hit families (Step 34d, best E <= 1e-3 vs PDB or CATH50) the hit names and a
#   proposed concordance call (for manual review).
# geNomad (results/step34j/mge/genomad, default thresholds): contig is "virus" if in contigs_virus_summary.tsv as a
#   whole contig, "plasmid" if in contigs_plasmid_summary.tsv; a member is in a provirus if its CDS lies within a
#   provirus region (virus summary rows named <contig>|provirus_<start>_<end>). Scores from the aggregated
#   classification table. DefenseFinder (results/step34j/mge/defensefinder): member protein in a system, or a system
#   protein within +-5 genes on the same contig ("defence"). Fractions are over distinct member contigs (members for
#   provirus and defence).
# Outputs (refuse to overwrite): results/step34j/context_all_families.tsv (4,509 families),
#   results/step34j_context.tsv (186 structural-hit families, Part A).
import os, sys, re, glob, collections
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts"); import step34j_themes as TH
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34j/"
OA, OC = D + "context_all_families.tsv", R + "step34j_context.tsv"
for p in (OA, OC):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
sc95 = {r["family"]: r for r in tsv(D + "context_scores.tsv") if r["set"] == "EU"}
sc = {r["family"]: r for r in tsv(D + "context_scores_nr50c.tsv") if r["set"] == "EU"}
jac = {r["family"]: r for r in tsv(D + "context_jaccard.tsv") if r["set"] == "EU"}
def jp(f):
    x = jac.get(f, {}).get("jaccard_percentile_vs_controls", "NA")
    return float(x) if x != "NA" else -1.0
for f, r in sc.items():
    r["conserved_context"] = str(int(int(r["n_nr_contexts"]) >= 3 and jp(f) >= 90))
ann = {r["rep95"]: r for r in tsv(D + "rep_annot.tsv")}
def atext(r95):
    a = ann.get(r95)
    return " ".join(a[k] for k in ("pfam", "kegg_ko", "preferred_name", "description") if a[k] != "-") if a else ""
seeds = collections.defaultdict(list)
for r in tsv(D + "context_nr50c_seeds.tsv"):
    if r["set"] == "EU": seeds[r["family"]].append(r["seed_member"])
nbt = collections.defaultdict(set)
for r in tsv(D + "context_neighbours.tsv"):
    if r["set"] == "EU": nbt[r["member"]] |= TH.themes(atext(r["rep95"]))
tfrac = {}
for f, ss in seeds.items():
    c = collections.Counter(t for s in ss for t in nbt.get(s, ()))
    tfrac[f] = {t: n / len(ss) for t, n in c.items()}
mem = collections.defaultdict(list)
for r in tsv(D + "context_members.tsv"):
    if r["set"] == "EU": mem[r["family"]].append(r)
G = D + "mge/genomad/"
def one(pat):
    x = glob.glob(G + pat)
    if len(x) != 1: raise SystemExit("ABORT: expected one file for %s, found %s" % (pat, x))
    return x[0]
vir, prov, pla, gsc = set(), collections.defaultdict(list), set(), {}
for r in tsv(one("*_summary/*_virus_summary.tsv")):
    n = r["seq_name"]
    if "|provirus_" in n:
        c, p = n.split("|provirus_"); a, b = p.split("_"); prov[c].append((int(a), int(b)))
    else: vir.add(n)
for r in tsv(one("*_summary/*_plasmid_summary.tsv")): pla.add(r["seq_name"])
for r in tsv(one("contigs_aggregated_classification/contigs_aggregated_classification.tsv")):
    gsc[r["seq_name"]] = (float(r["chromosome_score"]), float(r["plasmid_score"]), float(r["virus_score"]))
print("geNomad: virus contigs %d, contigs with provirus %d, plasmid contigs %d, scored %d" % (len(vir), len(prov), len(pla), len(gsc)))
dsys = glob.glob(D + "mge/defensefinder/*defense_finder_systems.tsv")
if len(dsys) != 1: raise SystemExit("ABORT: DefenseFinder systems table %s" % dsys)
insys = {}; near = collections.defaultdict(list)
for r in tsv(dsys[0]):
    for p in r["protein_in_syst"].split(","):
        insys[p] = r["subtype"]; c, n = p.rsplit("_", 1); near[c].append((int(n), r["subtype"]))
print("DefenseFinder systems %d, proteins in systems %d" % (sum(1 for _ in tsv(dsys[0])), len(insys)))
def mge(f):
    ms = mem.get(f, [])
    if not ms: return {}, "NA"
    ctg = {"%s__%s" % (m["stem"], m["contig"]) for m in ms}
    pv = df = 0; dsub = collections.Counter()
    for m in ms:
        c = "%s__%s" % (m["stem"], m["contig"]); a, b = int(m["start"]), int(m["end"])
        if any(x <= a and b <= y for x, y in prov.get(c, [])): pv += 1
        g = "%s_%s" % (c, m["gene_n"]); hit = insys.get(g) or next((s for n, s in near.get(c, []) if abs(n - int(m["gene_n"])) <= 5), None)
        if hit: df += 1; dsub[hit] += 1
    fr = {"virus": sum(c in vir for c in ctg) / len(ctg), "plasmid": sum(c in pla for c in ctg) / len(ctg),
          "provirus": pv / len(ms), "defence": df / len(ms)}
    vs = [gsc[c][2] for c in ctg if c in gsc]; ps = [gsc[c][1] for c in ctg if c in gsc]
    st = "contigs %d: virus %.2f plasmid %.2f | members %d: provirus %.2f defence %.2f%s | mean virus_score %s plasmid_score %s" % (
        len(ctg), fr["virus"], fr["plasmid"], len(ms), fr["provirus"], fr["defence"],
        (" (" + ",".join("%s:%d" % x for x in dsub.most_common(3)) + ")") if dsub else "",
        "%.3f" % (sum(vs) / len(vs)) if vs else "NA", "%.3f" % (sum(ps) / len(ps)) if ps else "NA")
    return fr, st
fams = [r["family_rep"] for r in tsv(R + "step34c2c_strict_v2_families.tsv") if r["recurrent"] == "1" and r["excluded"] == "0"]
pdbt = {}
for l in open(D + "ref_pdb/entries.idx", encoding="latin-1"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 4 and len(f[0]) == 4: pdbt[f[0].upper()] = f[3]
cdom, cname = {}, {}
for l in open(D + "ref_cath/cath-domain-list.txt"):
    if l[0] != "#": f = l.split(); cdom[f[0]] = ".".join(f[1:5])
for l in open(D + "ref_cath/cath-names.txt"):
    if l[0] != "#":
        f = l.rstrip("\n").split(None, 2)
        if len(f) == 3: cname[f[0]] = f[2].lstrip(":")
best = {}
for db in ("PDB", "CATH50"):
    for l in open(R + "step34d_prostt5/queries_vs_%s.tsv" % db):
        x = l.rstrip("\n").split("\t"); e = float(x[2])
        if e > 1e-3: continue
        k = (x[0], db)
        if k not in best or (e, -float(x[3])) < (best[k][1], -best[k][2]): best[k] = (x[1], e, float(x[3]))
qmap = {r["family_rep"]: r["query"] for r in tsv(R + "step34d_query_map.tsv")}
def hitname(db, t):
    if db == "PDB": return pdbt.get(t[:4].upper(), "NA")
    code = t.rsplit("_", 1)[1] if t.startswith("af_") else cdom.get(t, "NA")
    return "%s %s" % (code, cname.get(code, "NA"))
H = ["family_rep", "query", "n_strict_members", "n_members_with_context", "n_members_no_neighbour", "n_nr_contexts", "frac_nr_annotated",
     "top_neighbour50", "conservation", "percentile_vs_controls", "n_nr_contexts_rep95", "mean_pairwise_jaccard", "jaccard_percentile_vs_controls", "conserved_context", "top_neighbour_annotation", "strand_conservation",
     "distance_conservation", "median_gap_bp", "top_pfam", "top_kegg_ko", "top_kegg_module", "context_themes", "genomad_defensefinder",
     "frac_virus_contigs", "frac_plasmid_contigs", "frac_members_provirus", "frac_members_defence"]
HH = ["prostt5_hit", "best_pdb", "best_pdb_evalue", "best_pdb_title", "best_cath", "best_cath_evalue", "best_cath_superfamily",
      "concordance_proposed", "concordance_evidence"]
na, nc = 0, 0; calls = collections.Counter()
with open(OA, "w") as oa, open(OC, "w") as oc:
    oa.write("\t".join(H + HH) + "\n"); oc.write("\t".join(H + HH) + "\n")
    for f in fams:
        s = sc.get(f, {}); q = qmap[f]
        fr, st = mge(f)
        th = tfrac.get(f, {})
        row = [f, q] + [s.get(k, "NA") for k in H[2:6]] + [s.get("frac_nr_annotated", "NA"), s.get("top_neighbour50", "NA"), s.get("conservation", "NA"),
               s.get("percentile_vs_controls", "NA"), sc95.get(f, {}).get("n_nr_contexts", "NA"), jac.get(f, {}).get("mean_pairwise_jaccard", "NA"),
               jac.get(f, {}).get("jaccard_percentile_vs_controls", "NA"), s.get("conserved_context", "0"), s.get("top_neighbour_annotation", "NA"),
               s.get("strand_conservation", "NA"), s.get("distance_conservation", "NA"), s.get("median_gap_bp", "NA"),
               s.get("top_pfam", "NA"), s.get("top_kegg_ko", "NA"), s.get("top_kegg_module", "NA"),
               "|".join("%s:%.2f" % (t, v) for t, v in sorted(th.items(), key=lambda x: -x[1])) or "-", st] + \
              ["%.3f" % fr[k] if fr else "NA" for k in ("virus", "plasmid", "provirus", "defence")]
        bp, bc = best.get((q, "PDB")), best.get((q, "CATH50"))
        hit = int(bool(bp or bc))
        if hit:
            text = " ".join(x for x in ((hitname("PDB", bp[0]) if bp else ""), (hitname("CATH50", bc[0]) if bc else "")) if x)
            nnr = int(s.get("n_nr_contexts", 0) or 0)
            call, ev = TH.concordance(text, th, nnr, s.get("conserved_context") == "1", s.get("top_neighbour_annotation", "-"), fr)
            calls[call] += 1
        else: call, ev = "NA", "-"
        hh = [hit] + ([bp[0], "%.2e" % bp[1], hitname("PDB", bp[0])] if bp else ["-"] * 3) + ([bc[0], "%.2e" % bc[1], hitname("CATH50", bc[0])] if bc else ["-"] * 3) + \
             ["proposed: " + call if hit else "NA", ev]
        line = "\t".join(map(str, row + hh)).replace("\n", " ") + "\n"
        oa.write(line); na += 1
        if hit: oc.write(line); nc += 1
print("families %d; structural-hit families %d; proposed calls %s" % (na, nc, dict(calls)))
if nc != 186: raise SystemExit("ABORT: expected 186 structural-hit families, got %d" % nc)
print("wrote", OA, OC)
