# Step 34j D1: summaries of the v2 structure classification (scripts/step34j_b4_classify.py, 2026-10-23 run).
# Sets: EU = models of the 4,250 member-validated recurrent EU families; control = 4,509 length-matched known reps.
# pLDDT bin = mean C-alpha pLDDT of the ESMFold model: >= 70, 50-70, < 50.
# Hit = Foldseek E <= 1e-3 (results/step34j/foldseek/<db>.m8, full model or any of its Merizo domains); the stricter
#   criterion adds qTM >= 0.5. AFDB50_annotated = a hit among the 25 best AFDB50 hits of its query unit whose UniProt
#   accession is annotated (results/step34j/uniprot_annotation_v2.tsv), under both ProtNLM rules.
# Outputs (results/step34j/; refuse to overwrite):
#   models_by_set_plddt_v2.tsv   models per set and pLDDT bin, with hit class counts (both ProtNLM rules)
#   hit_share_by_db_v2.tsv       per database and any database, models with a hit, by set and pLDDT bin, both criteria
#   confident_eu_no_hit_v2.tsv   EU models with pLDDT >= 70 and no hit (E <= 1e-3) in any of the four databases
#   merizo_v2.tsv                Merizo result per model >= 200 aa, with hits of the full model vs domains only
#   merizo_summary_v2.tsv        domain counts and domain-only hits by set
#   search_settings.tsv          per database: version, entries, Foldseek command settings, coverage, hit thresholds
import os, re, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34j/"
OUT = {k: D + k + ".tsv" for k in ("models_by_set_plddt_v2", "hit_share_by_db_v2", "confident_eu_no_hit_v2", "merizo_v2", "merizo_summary_v2", "search_settings")}
for p in OUT.values():
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
S = {r["id"]: r for r in tsv(R + "step34j_structures_v2.tsv")}
ann = {r["accession"]: r for r in tsv(D + "uniprot_annotation_v2.tsv")}
SETS = ("EU", "control"); PB = (">=70", "50-70", "<50")
n_set = collections.Counter(r["set"] for r in S.values()); print("models per set", dict(n_set), "no_model", sum(r["hit_class"] == "no_model" for r in S.values()))
units = set(S) | {f[:-4] for f in os.listdir(D + "domains")}
def unit(q):
    for c in (q, re.sub(r"\.pdb$", "", q), re.sub(r"(\.pdb)?_[A-Z]$", "", q)):
        if c in units: return c
    raise SystemExit("ABORT: cannot map query %s" % q)
def model_of(u): return u.split("__D")[0] if re.search(r"__D\d+$", u) else u
def acc(t):
    m = re.match(r"AF-([A-Z0-9]+)-F\d+", t); return m.group(1) if m else None
DBS = ("PDB", "CATH50", "AFDB50", "ESMAtlas30")
H = {db: collections.defaultdict(list) for db in DBS}
for db in DBS:
    for l in open(D + "foldseek/%s.m8" % db):
        f = l.rstrip("\n").split("\t"); e = float(f[6])
        if e > 1e-3: continue
        u = unit(f[0]); H[db][model_of(u)].append(dict(unit=u, t=f[1], e=e, bits=float(f[7]), qtm=float(f[10])))
for hs in H["AFDB50"].values():
    byu = collections.defaultdict(list)
    for h in hs: byu[h["unit"]].append(h)
    for v in byu.values():
        for h in sorted(v, key=lambda h: (h["e"], -h["bits"]))[:25]: h["top25"] = True
def afa(h, rule): a = acc(h["t"]); return bool(h.get("top25")) and a in ann and ann[a]["annotated_protnlm_" + rule] == "1"
def has(i, db, crit, rule="unannotated"):
    q = (lambda h: h["qtm"] >= 0.5) if crit == "E<=1e-3,qTM>=0.5" else (lambda h: True)
    if db == "experimental": return any(q(h) for d in ("PDB", "CATH50") for h in H[d].get(i, []))
    if db == "AFDB50_annotated": return any(q(h) and afa(h, rule) for h in H["AFDB50"].get(i, []))
    if db == "any": return any(q(h) for d in DBS for h in H[d].get(i, []))
    return any(q(h) for h in H[db].get(i, []))
CL = ["experimental", "AFDB_annotated", "AFDB_unannotated_only", "ESMAtlas_only", "none", "no_model"]
with open(OUT["models_by_set_plddt_v2"], "w") as o:
    o.write("protnlm_rule\tset\tplddt_bin\tn_models\t" + "\t".join(CL) + "\tstructure_confident_annotated\n")
    for rule, cc, fc in (("unannotated", "hit_class", "structure_confident_annotated"), ("annotated", "hit_class_protnlm_annotated", "structure_confident_protnlm_annotated")):
        for s in SETS:
            for b in PB + ("all",):
                X = [r for r in S.values() if r["set"] == s and (b == "all" or r["plddt_bin"] == b)]; c = collections.Counter(r[cc] for r in X)
                o.write("%s\t%s\t%s\t%d\t%s\t%d\n" % (rule, s, b, len(X), "\t".join(str(c[k]) for k in CL), sum(r[fc] == "1" for r in X)))
DBR = ("PDB", "CATH50", "experimental", "AFDB50", "AFDB50_annotated", "ESMAtlas30", "any")
with open(OUT["hit_share_by_db_v2"], "w") as o:
    o.write("criterion\tdatabase\tset\tplddt_bin\tn_models\tn_with_hit\tpct\n")
    for crit in ("E<=1e-3", "E<=1e-3,qTM>=0.5"):
        for db in DBR:
            for s in SETS:
                for b in PB + ("all",):
                    X = [i for i, r in S.items() if r["set"] == s and r["hit_class"] != "no_model" and (b == "all" or r["plddt_bin"] == b)]
                    k = sum(has(i, db, crit) for i in X)
                    o.write("%s\t%s\t%s\t%s\t%d\t%d\t%.2f\n" % (crit, db, s, b, len(X), k, 100.0 * k / len(X) if X else 0))
                    if b == "all": print("%-17s %-17s %-7s %4d / %4d (%.2f%%)" % (crit, db, s, k, len(X), 100.0 * k / len(X)))
    # AFDB50_annotated under the ProtNLM-accepted rule, for comparison
    for crit in ("E<=1e-3", "E<=1e-3,qTM>=0.5"):
        for s in SETS:
            for b in PB + ("all",):
                X = [i for i, r in S.items() if r["set"] == s and r["hit_class"] != "no_model" and (b == "all" or r["plddt_bin"] == b)]
                k = sum(has(i, "AFDB50_annotated", crit, "accepted") for i in X)
                o.write("%s\tAFDB50_annotated_protnlm_accepted\t%s\t%s\t%d\t%d\t%.2f\n" % (crit, s, b, len(X), k, 100.0 * k / len(X) if X else 0))
NH = [i for i, r in sorted(S.items()) if r["set"] == "EU" and r["plddt_bin"] == ">=70" and not has(i, "any", "E<=1e-3")]
with open(OUT["confident_eu_no_hit_v2"], "w") as o:
    o.write("id\tfamily_rep\tlength_aa\tmean_ca_plddt\tptm\tn_merizo_domains\n")
    for i in NH: r = S[i]; o.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (i, r["family_rep"], r["length_aa"], r["mean_ca_plddt"], r["ptm"], r["n_merizo_domains"]))
print("EU models pLDDT >= 70: %d; with no hit (E <= 1e-3) in any database: %d" % (sum(r["set"] == "EU" and r["plddt_bin"] == ">=70" for r in S.values()), len(NH)))
mz = {}
for l in open(D + "merizo/merizo_results.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "input": continue
    mz[f[0][:-4]] = f
ge = [l.strip() for l in open(D + "merizo/ids_ge200.txt") if l.strip()]
ge = [x[:-4] if x.endswith(".pdb") else x for x in ge]
nomod = [x for x in ge if S[x]["hit_class"] == "no_model"]; miss = [x for x in ge if x not in mz and x not in nomod]
print("Merizo: inputs >= 200 aa %d, without an ESMFold model %d, results %d, missing %d" % (len(ge), len(nomod), len(mz), len(miss)))
if miss: raise SystemExit("ABORT: Merizo results missing for %d models" % len(miss))
def full_dom(i):
    fh = any(h["unit"] == i for db in DBS for h in H[db].get(i, [])); dh = any(h["unit"] != i for db in DBS for h in H[db].get(i, []))
    return fh, dh
with open(OUT["merizo_v2"], "w") as o:
    o.write("id\tset\tlength_aa\tmean_ca_plddt\tn_domains\tpIoU\tresidues_in_domains\tdomain_ranges\thit_full_model\thit_domain\thit_domain_only\thit_class\n")
    for i in sorted(mz):
        f = mz[i]; r = S[i]; fh, dh = full_dom(i)
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\n" % (i, r["set"], r["length_aa"], r["mean_ca_plddt"], f[4], f[5], f[2], f[7], fh, dh, dh and not fh, r["hit_class"]))
with open(OUT["merizo_summary_v2"], "w") as o:
    o.write("set\tn_models_ge200\tn_domains_0\tn_domains_1\tn_domains_2\tn_domains_3plus\ttotal_domains\tmodels_hit_full\tmodels_hit_domain\tmodels_hit_domain_only\n")
    for s in SETS + ("all",):
        X = [i for i in mz if s == "all" or S[i]["set"] == s]; nd = collections.Counter(min(int(mz[i][4]), 3) for i in X); fd = [full_dom(i) for i in X]
        row = (s, len(X), nd[0], nd[1], nd[2], nd[3], sum(int(mz[i][4]) for i in X), sum(a for a, _ in fd), sum(b for _, b in fd), sum(b and not a for a, b in fd))
        o.write("\t".join(map(str, row)) + "\n"); print("Merizo %-7s models %d, 0/1/2/3+ domains %d/%d/%d/%d, domains %d; hit full %d, domain %d, domain only %d" % row)
V = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(D + "foldseek_db_versions.tsv")}
DV = "/bigdata/stajichlab/lshad003/software/foldseek_dbs/"
pdbv = open(DV + "PDB.version").read().split()
with open(OUT["search_settings"], "w") as o:
    o.write("database\tversion\tn_entries\tfoldseek\tsearch_mode\tsearch_evalue\tcoverage\thit_rule\tconfident_rule\tqueries\n")
    common = ("foldseek commit 463739e0014a1549a527de589102cde98f802f37 (search -a; convertalis query,target,fident,alnlen,qlen,tlen,evalue,bits,prob,alntmscore,qtmscore,ttmscore,lddt,qcov,tcov)",
              "-e 10", "no coverage filter (-c 0, --cov-mode 0, Foldseek defaults; qcov kept)", "E <= 1e-3", "qTM >= 0.5, prob >= 0.5, model mean C-alpha pLDDT >= 70",
              "8,759 ESMFold models + Merizo domains of models >= 200 aa")
    rows = [("PDB", "PDB100 Foldseek build, PDB_DATE %s (FOLDSEEK_COMMIT %s)" % (pdbv[2], pdbv[4]), "NA", "--exhaustive-search 1"),
            ("CATH50", "CATH50 %s (Foldseek database); superfamily names CATH v4.4.0" % open(DV + "CATH50.version").read().strip(), "NA", "--exhaustive-search 1"),
            ("AFDB50", "Alphafold/UniProt50 %s" % " ".join(V["AFDB50"][2].split()), V["AFDB50"][5], "-s 9.5 --max-seqs 1000"),
            ("ESMAtlas30", "ESMAtlas30 %s" % " ".join(V["ESMAtlas30"][2].split()), V["ESMAtlas30"][5], "-s 9.5 --max-seqs 1000")]
    for db, v, n, mode in rows:
        o.write("\t".join([db, v, n, common[0], mode, common[1], common[2], common[3] + (" (AFDB annotation: 25 best hits per query unit; UniProt live, secondary or UniSave last version; ProtNLM-only names unannotated)" if db == "AFDB50" else ""), common[4], common[5]]) + "\n")
print("wrote", " ".join(OUT.values()))
