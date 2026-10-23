# Step 34j B4 UniProt gap: b4 (job 30171185) asked UniProt for 88,364 AFDB accessions and got 31,660 back.
# Rebuilds the requested set exactly as scripts/step34j_b4_classify.py (E <= 1e-3, 25 best AFDB50 hits per query unit),
#   reads the saved responses (results/step34j/uniprot/batch_*.tsv), and measures what the missing ones can change in
#   results/step34j_structures.tsv: a missing accession is treated as unannotated by b4, so only models whose hit_class
#   is AFDB_unannotated_only (no PDB/CATH50 hit) can move up; the upper bound assumes every missing accession is annotated.
# No network access. Outputs in results/step34j/uniprot_gap/: missing_accessions.tsv, per_batch.tsv, impact.tsv.
import os, re, glob, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; D = W + "/results/step34j/"; O = D + "uniprot_gap/"
os.makedirs(O, exist_ok=True)
units = {os.path.basename(x)[:-4] for x in glob.glob(D + "esmfold_pdb/*.pdb")} | {os.path.basename(x)[:-4] for x in glob.glob(D + "domains/*.pdb")}
def unit(q):
    for c in (q, re.sub(r"\.pdb$", "", q), re.sub(r"(\.pdb)?_[A-Z]$", "", q)):
        if c in units: return c
    raise SystemExit("ABORT: cannot map query %s" % q)
def model_of(u): return u.split("__D")[0] if re.search(r"__D\d+$", u) else u
def acc(t):
    m = re.match(r"AF-([A-Z0-9]+)-F\d+", t); return m.group(1) if m else None
byu = collections.defaultdict(list)
for l in open(D + "foldseek/AFDB50.m8"):
    f = l.rstrip("\n").split("\t"); e = float(f[6])
    if e > 1e-3: continue
    byu[unit(f[0])].append((e, -float(f[7]), f[1], float(f[10]), float(f[8])))
top = collections.defaultdict(list)   # model -> top-25 hits (acc, qtm, prob)
for u, v in byu.items():
    for e, nb, t, qtm, prob in sorted(v)[:25]: top[model_of(u)].append((acc(t), qtm, prob))
want = sorted({a for v in top.values() for a, _, _ in v if a})
print("requested accessions %d" % len(want))
got = set(); per = []
for bi in range(0, len(want), 500):
    f = D + "uniprot/batch_%05d.tsv" % (bi // 500)
    r = [l.split("\t")[0] for l in list(open(f))[1:]]
    req = set(want[bi:bi + 500]); got |= set(r)
    per.append((bi // 500, len(req), len(r), len(set(r) - req), want[bi], want[min(bi + 499, len(want) - 1)]))
miss = [a for a in want if a not in got]; extra = got - set(want)
print("returned %d; missing %d; returned but not requested (secondary->primary remaps) %d" % (len(got & set(want)), len(miss), len(extra)))
with open(O + "per_batch.tsv", "w") as o:
    o.write("batch\trequested\treturned\treturned_not_requested\tfirst_acc\tlast_acc\n")
    for p in per: o.write("\t".join(map(str, p)) + "\n")
fmt = collections.Counter("A0A-prefixed TrEMBL (10 char)" if a.startswith("A0A") else "6 char" if len(a) == 6 else "other %d char" % len(a) for a in miss)
fmt_all = collections.Counter("A0A-prefixed TrEMBL (10 char)" if a.startswith("A0A") else "6 char" if len(a) == 6 else "other %d char" % len(a) for a in want)
print("missing by format:", dict(fmt)); print("requested by format:", dict(fmt_all))
with open(O + "missing_accessions.tsv", "w") as o:
    o.write("accession\n"); o.writelines(a + "\n" for a in miss)
M = set(miss)
S = {}
with open(W + "/results/step34j_structures.tsv") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
    for l in fh:
        f = l.rstrip("\n").split("\t"); S[f[0]] = f
rows = []
for s in ("EU", "control"):
    c = collections.Counter(); 
    for i, f in S.items():
        if f[ix["set"]] != s: continue
        cls = f[ix["hit_class"]]; v = top.get(i, [])
        nm = sum(1 for a, _, _ in v if a in M)
        if v: c["models_with_AFDB_top25"] += 1
        if nm: c["models_with_a_missing_acc"] += 1
        if nm and len(v) == nm: c["models_all_top25_missing"] += 1
        if cls == "AFDB_annotated" and nm: c["AFDB_annotated_with_missing"] += 1
        if cls == "AFDB_unannotated_only":
            c["AFDB_unannotated_only"] += 1
            if nm:
                c["unannot_only_could_move_to_AFDB_annotated"] += 1
                p = f[ix["mean_ca_plddt"]]
                if p != "NA" and float(p) >= 70 and any(a in M and q >= 0.5 and pr >= 0.5 for a, q, pr in v):
                    c["could_become_structure_confident"] += 1
        if cls == "AFDB_annotated" and f[ix["best_db"]] == "AFDB50" and f[ix["best_target_name"]] == "NA": c["AFDB_annotated_best_name_NA"] += 1
        if cls == "AFDB_unannotated_only" and f[ix["best_target_name"]] == "NA": c["unannot_only_best_name_NA"] += 1
        c["structure_confident_now"] += int(f[ix["structure_confident_annotated"]])
    for k in ("models_with_AFDB_top25", "models_with_a_missing_acc", "models_all_top25_missing", "AFDB_annotated_with_missing",
              "AFDB_unannotated_only", "unannot_only_could_move_to_AFDB_annotated", "unannot_only_best_name_NA",
              "structure_confident_now", "could_become_structure_confident"):
        rows.append((s, k, c[k])); print("%-7s %-45s %d" % (s, k, c[k]))
with open(O + "impact.tsv", "w") as o:
    o.write("set\tmetric\tn\n"); o.writelines("%s\t%s\t%d\n" % r for r in rows)
print("DONE")
