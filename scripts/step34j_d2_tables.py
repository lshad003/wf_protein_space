# Step 34j D2: supplementary Tables S5, S6 and S7.
# Table S5 (results/supp_table_s5_structures.tsv): one row per Step 34j ESMFold input (4,250 EU family models, 4,509
#   controls; one EU input has no model), from results/step34j_structures_v2.tsv. hit_class and
#   structure_confident_annotated count ProtNLM-only names (ECO:0008006) as unannotated; the *_protnlm_annotated columns
#   accept them.
# Table S6 (results/supp_fold_candidates.tsv; the Step 29 version, built from 95% cluster members, is kept as
#   results/supp_fold_candidates_step29.tsv): the representatives of the three Step 29 candidate families
#   (UHM626.41025__k141_276632_94, UHM634.41029__k141_110414_38, UHM22.35747__k141_101893_2; Step 29 member models
#   UHM125.35769__k141_151130_2, UHM145.35772__k141_106404_54, UHM145.35771__k141_106979_4, results/clusters_95_n88.tsv)
#   with their Step 34j models and best hit per database (E <= 1e-3), and Cse2 variant A (UHM45.35758__k141_100597_236;
#   its family is class_any K in results/step34h_member_classes.tsv, so it has no Step 34j model; Step 34h model and
#   PDB/CATH50 hits from results/step34h_cse2_varA.tsv, no AFDB50/ESMAtlas30 search). Context columns from
#   results/step34j_matrix_v2.tsv (family level).
# Table S7 (results/supp_table_s7_c_id95.tsv): results/step34s/side_by_side.tsv with c_id95 (primary) first; the
#   unfiltered values are kept as superseded.
# Refuses to overwrite S5, S7 and supp_fold_candidates_step29.tsv; S6 is replaced only after the Step 29 copy exists.
import os, re, shutil, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34j/"
S5, S6, S6OLD, S7 = R + "supp_table_s5_structures.tsv", R + "supp_fold_candidates.tsv", R + "supp_fold_candidates_step29.tsv", R + "supp_table_s7_c_id95.tsv"
for p in (S5, S7, S6OLD):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
S = {r["id"]: r for r in tsv(R + "step34j_structures_v2.tsv")}
# S5
C5 = ["id", "set", "family_rep", "length_aa", "mean_ca_plddt", "esmfold_plddt", "ptm", "plddt_bin", "n_merizo_domains", "hit_class",
      "hit_class_protnlm_annotated", "structure_confident_annotated", "structure_confident_protnlm_annotated", "best_db", "best_target",
      "best_evalue", "best_target_name", "confident_db", "confident_target", "confident_qtm", "confident_target_name", "prostt5_hit"]
for db in ("PDB", "CATH50", "AFDB50", "ESMAtlas30"): C5 += ["%s_%s" % (db, x) for x in ("target", "evalue", "qtm", "prob")]
C5 += ["afdb_annot_sources"]
with open(S5, "w") as o:
    o.write("\t".join(C5) + "\n")
    for i in sorted(S, key=lambda i: (S[i]["set"] != "EU", i)): o.write("\t".join(S[i][c] for c in C5) + "\n")
print("Table S5: %d rows (%s)" % (len(S), dict(collections.Counter(r["set"] for r in S.values()))))
# names for per-database best hits
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
ua = {r["accession"]: r for r in tsv(D + "uniprot_annotation_v2.tsv")}
def tname(db, t):
    if t in ("-", ""): return "-"
    if db == "PDB": return pdbt.get(t[:4].upper(), "NA")
    if db == "CATH50":
        code = t.rsplit("_", 1)[1] if t.startswith("af_") else cdom.get(t, "NA"); return "%s %s" % (code, cname.get(code, "NA"))
    if db == "AFDB50":
        m = re.match(r"AF-([A-Z0-9]+)-F\d+", t); a = ua.get(m.group(1)) if m else None
        return "NA" if a is None else "%s [%s; annotated %s]" % (a["name"], a["source"], a["annotated_protnlm_unannotated"])
    return "-"
# S6
if not os.path.exists(S6): raise SystemExit("ABORT: " + S6 + " missing")
shutil.copy2(S6, S6OLD)
old = {r["rep"]: r for r in tsv(S6OLD)}
MX = {r["family_rep"]: r for r in tsv(R + "step34j_matrix_v2.tsv")}
CAND = [("UHM45.35758__k141_100597_236", "UHM636.41030__k141_332084_5"), ("UHM626.41025__k141_276632_94", "UHM125.35769__k141_151130_2"),
        ("UHM634.41029__k141_110414_38", "UHM145.35772__k141_106404_54"), ("UHM22.35747__k141_101893_2", "UHM145.35771__k141_106979_4")]
H6 = ["rep", "model_source", "step29_model", "step29_label", "length_aa", "mean_ca_plddt", "ptm", "hit_class", "hit_class_protnlm_annotated",
      "structure_confident_annotated", "structure_confident_protnlm_annotated"]
for db in ("PDB", "CATH50", "AFDB50", "ESMAtlas30"): H6 += ["%s_best_target" % db, "%s_evalue" % db, "%s_qtm" % db, "%s_name" % db]
H6 += ["n_nr_contexts", "conserved_context", "top_neighbour_annotation", "c1_rank", "c1_rank_tier"]
with open(S6, "w") as o:
    o.write("\t".join(H6) + "\n")
    for rep, m29 in CAND:
        x = MX[rep]; ctx = [x["n_nr_contexts"], x["conserved_context"], x["top_neighbour_annotation"], x["rank"], x["rank_tier"]]
        if rep in S:
            r = S[rep]; row = [rep, "Step 34j ESMFold (results/step34j/esmfold_pdb)", m29, old[m29]["candidate"], r["length_aa"], r["mean_ca_plddt"], r["ptm"],
                              r["hit_class"], r["hit_class_protnlm_annotated"], r["structure_confident_annotated"], r["structure_confident_protnlm_annotated"]]
            for db in ("PDB", "CATH50", "AFDB50", "ESMAtlas30"):
                t = r["%s_target" % db]; row += [t, r["%s_evalue" % db], r["%s_qtm" % db], tname(db, t)]
        else:
            h = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(R + "step34h_cse2_varA.tsv")}
            def top(k):
                f = h[k].split(); return [f[0], "%.2e" % float(f[2]), "%.3f" % float(f[6])]
            p, c = top("varA PDB top1"), top("varA CATH50 top1")
            row = [rep, "Step 34h ESMFold (results/esmfold_pdb); not in Step 34j (family class_any K)", m29, old[m29]["candidate"], h["residues"],
                   h["mean pLDDT, CA atoms"], h["esm-fold reported pLDDT / pTM"].split("/")[1].strip(), "NA", "NA", "NA", "NA"]
            row += p + [tname("PDB", p[0])] + c + [tname("CATH50", re.sub(r"^(\w{4})(\w)(\d\d)$", r"\1\2\3", c[0]))] + ["not searched"] * 8
        o.write("\t".join(str(v).replace("\t", " ") for v in row + ctx) + "\n")
        print("Table S6 %s: %s" % (rep, " | ".join("%s=%s" % (k, v) for k, v in zip(H6[4:], (row + ctx)[4:]))))
# S7
with open(S7, "w") as o:
    o.write("section\titem\tc_id95\tc_id90\tc_id95_q1\tunfiltered_superseded\tsource_file_in_results_step34s_label\n")
    n = 0
    for r in tsv(R + "step34s/side_by_side.tsv"):
        o.write("\t".join([r["section"], r["item"], r["c_id95"], r["c_id90"], r["c_id95_q1"], r["old"], r["file (in results/step34s/<label>/)"]]) + "\n"); n += 1
print("Table S7: %d rows" % n)
print("wrote", S5, S6, S6OLD, S7)
