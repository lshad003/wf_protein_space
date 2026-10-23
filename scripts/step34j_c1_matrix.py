# Step 34j C1: evidence matrix (one row per recurrent EU family, 4,509) and proposed candidate ranking (top 20).
# Inputs: results/step34c2c_strict_v2_families.tsv (n_variants95, n_member_proteins), results/step34h_member_classes.tsv
#   (member-validated = every 95% cluster class_any == EU; strict unique variants; rep_partial of the query rep),
#   results/step34d_query_map.tsv, results/step34j/context_all_families.tsv (Part A), results/step34j_structures.tsv
#   (Part B; member-validated families only, others NA).
# Context columns: n_nr_contexts = A2c backgrounds; conserved_context as step34j_a4 (>= 3 backgrounds, Jaccard pct >= 90).
# Concordance for the coordinate hit: scripts/step34j_themes.py on the confident hit name (else the best hit name of an
#   experimental or AFDB-annotated class). Ranking tiers (proposed, no function assigned):
#   1 structure-confident annotated hit + conserved context + concordance supports (coordinate or ProstT5 hit)
#   2 structure-confident annotated hit + supports, context not conserved
#   3 structure-confident annotated hit + conserved context, neutral
#   4 structure-confident annotated hit only
#   5 any experimental or AFDB-annotated hit, or a ProstT5 hit
#   6 the rest
#   Within a tier: Jaccard context percentile vs controls (desc), NR contexts (desc), confident qTM (desc), family ID.
# Outputs (refuse to overwrite): results/step34j_matrix.tsv, results/step34j_top20.tsv
# 2026-10-22: reads results/step34j_structures_v2.tsv (b4 with the UniProt gap closed) and writes step34j_matrix_v2.tsv,
#   step34j_top20_v2.tsv; run 1 (job 30171186) outputs step34j_matrix.tsv, step34j_top20.tsv kept.
import os, sys, collections
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts"); import step34j_themes as TH
R = "/bigdata/stajichlab/lshad003/wf_protein_space/results/"; OM, OT = R + "step34j_matrix_v2.tsv", R + "step34j_top20_v2.tsv"
for p in (OM, OT):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
fam = {r["family_rep"]: r for r in tsv(R + "step34c2c_strict_v2_families.tsv") if r["recurrent"] == "1" and r["excluded"] == "0"}
mc = collections.defaultdict(list); part = {}
for r in tsv(R + "step34h_member_classes.tsv"):
    if r["family_rep"] in fam: mc[r["family_rep"]].append(r); part[r["rep"]] = r["rep_partial"]
qm = {r["family_rep"]: r for r in tsv(R + "step34d_query_map.tsv")}
ctx = {r["family_rep"]: r for r in tsv(R + "step34j/context_all_families.tsv")}
st = {r["family_rep"]: r for r in tsv(R + "step34j_structures_v2.tsv") if r["set"] == "EU"}
if len(fam) != 4509 or len(ctx) != 4509: raise SystemExit("ABORT: family counts %d %d" % (len(fam), len(ctx)))
def fl(x, d=-1.0):
    try: return float(x)
    except (TypeError, ValueError): return d
rows = []
for f in sorted(fam):
    c, s, q = ctx[f], st.get(f, {}), qm[f]
    mv = int(all(r["class_any"] == "EU" for r in mc[f]))
    nuv = sum(int(r["strict_unique_variants"]) for r in mc[f])
    cls = s.get("hit_class", "not_modelled")
    conf = s.get("structure_confident_annotated", "0") == "1"
    htxt = s.get("confident_target_name") if conf else (s.get("best_target_name") if cls in ("experimental", "AFDB_annotated") else None)
    th = {}
    for x in (c["context_themes"].split("|") if c["context_themes"] != "-" else []):
        t, v = x.rsplit(":", 1); th[t] = float(v)
    fr = {k: fl(c["frac_" + n], 0) for k, n in (("virus", "virus_contigs"), ("plasmid", "plasmid_contigs"), ("provirus", "members_provirus"), ("defence", "members_defence"))}
    if htxt and htxt not in ("-", "NA"):
        cc, ce = TH.concordance(htxt, th, int(fl(c["n_nr_contexts"], 0)), c["conserved_context"] == "1", c["top_neighbour_annotation"], fr)
    else: cc, ce = "NA", "-"
    cons = c["conserved_context"] == "1"
    sup = cc == "supports" or c["concordance_proposed"] == "proposed: supports"
    anyhit = cls in ("experimental", "AFDB_annotated") or c["prostt5_hit"] == "1"
    tier = 1 if conf and cons and sup else 2 if conf and sup else 3 if conf and cons else 4 if conf else 5 if anyhit else 6
    rows.append(dict(family_rep=f, query=q["query"], member_validated=mv, length_aa=q["length_aa"], query_partial=part.get(q["query"], "NA"),
        n_95_variants=fam[f]["n_variants95"], n_strict_unique_variants=nuv, n_member_proteins=fam[f]["n_member_proteins"],
        n_metagenomes=fam[f]["n_metagenomes"], n_nr_contexts=c["n_nr_contexts"], prostt5_hit=c["prostt5_hit"],
        prostt5_best=(c["best_pdb_title"] if c["best_pdb"] != "-" else c["best_cath_superfamily"]),
        hit_class=cls, best_db=s.get("best_db", "NA"), best_target=s.get("best_target", "NA"), best_evalue=s.get("best_evalue", "NA"),
        best_target_name=s.get("best_target_name", "NA"), structure_confident_annotated=int(conf), confident_db=s.get("confident_db", "NA"),
        confident_target=s.get("confident_target", "NA"), confident_qtm=s.get("confident_qtm", "NA"), confident_target_name=s.get("confident_target_name", "NA"),
        mean_ca_plddt=s.get("mean_ca_plddt", "NA"), ptm=s.get("ptm", "NA"), n_merizo_domains=s.get("n_merizo_domains", "NA"),
        n_nr_contexts_rep95=c["n_nr_contexts_rep95"], conservation=c["conservation"], percentile_vs_controls=c["percentile_vs_controls"],
        mean_pairwise_jaccard=c["mean_pairwise_jaccard"], jaccard_percentile_vs_controls=c["jaccard_percentile_vs_controls"], conserved_context=c["conserved_context"],
        top_neighbour_annotation=c["top_neighbour_annotation"], strand_conservation=c["strand_conservation"], distance_conservation=c["distance_conservation"],
        top_pfam=c["top_pfam"], top_kegg_ko=c["top_kegg_ko"], top_kegg_module=c["top_kegg_module"], context_themes=c["context_themes"],
        genomad_defensefinder=c["genomad_defensefinder"], concordance_prostt5=c["concordance_proposed"], concordance_prostt5_evidence=c["concordance_evidence"],
        concordance_coordinate=("proposed: " + cc) if cc != "NA" else "NA", concordance_coordinate_evidence=ce, rank_tier=tier))
rows.sort(key=lambda r: (r["rank_tier"], -fl(r["jaccard_percentile_vs_controls"]), -fl(r["n_nr_contexts"]), -fl(r["confident_qtm"]), r["family_rep"]))
for i, r in enumerate(rows, 1): r["rank"] = i
H = list(rows[0].keys())
with open(OM, "w") as o:
    o.write("\t".join(H) + "\n")
    for r in rows: o.write("\t".join(str(r[k]).replace("\t", " ").replace("\n", " ") for k in H) + "\n")
def summ(r):
    hit = ("%s %s (%s; qTM %s)" % (r["confident_db"], r["confident_target"], r["confident_target_name"][:70], r["confident_qtm"])
           if r["structure_confident_annotated"] else "%s %s (%s)" % (r["best_db"], r["best_target"], r["best_target_name"][:70]))
    return "%s aa, pLDDT %s, pTM %s; %s; context %s backgrounds, Jaccard %s (pct %s), top neighbour %s in %s; coord %s; ProstT5 %s; %s" % (
        r["length_aa"], r["mean_ca_plddt"], r["ptm"], hit, r["n_nr_contexts"], r["mean_pairwise_jaccard"], r["jaccard_percentile_vs_controls"],
        r["top_neighbour_annotation"][:60], r["conservation"], r["concordance_coordinate"], r["concordance_prostt5"], r["genomad_defensefinder"][:90])
with open(OT, "w") as o:
    o.write("rank\trank_tier\tfamily_rep\tquery\thit_class\tconcordance_coordinate\tconcordance_prostt5\tevidence_summary\n")
    for r in rows[:20]:
        o.write("%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (r["rank"], r["rank_tier"], r["family_rep"], r["query"], r["hit_class"], r["concordance_coordinate"], r["concordance_prostt5"], summ(r).replace("\t", " ")))
print("tiers:", dict(sorted(collections.Counter(r["rank_tier"] for r in rows).items())))
print("member-validated %d; structure-confident annotated %d; conserved context %d" % (sum(r["member_validated"] for r in rows),
      sum(r["structure_confident_annotated"] for r in rows), sum(r["conserved_context"] == "1" for r in rows)))
print("wrote", OM, OT)
