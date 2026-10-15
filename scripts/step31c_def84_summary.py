# Step 31c summary: per model of the 84 defence-context single-variant reps (results/def84_labels.tsv).
# pLDDT and pTM from the ESMFold logs: the latest logs/step31b_esmfold.*.log for the 73 new models and
#   logs/step29b_esmfold.29932169.log for the 11 Step 29 copies; plddt_ca (mean CA B-factor) is recomputed from every PDB
#   as a check. helix/strand fractions, Rg ratio and class exactly as scripts/step29_peek_v2.py (class: pLDDT >= 70 and
#   pTM >= 0.5 -> folded, or helix_only if helix >= 0.8 and strand < 0.05; otherwise low_confidence).
# Foldseek: top hit per database by bits (as step29e) from results/step31_foldseek_{PDB,CATH50}.m8.
# Output: results/step31_def84_summary.tsv
import glob, os, re
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; Q = R + "foldseek_def84"
OUT = R + "step31_def84_summary.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
lab = {}
with open(R + "def84_labels.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); lab[f["rep"]] = f
if len(lab) != 84: raise SystemExit("ABORT: expected 84 labels, got %d" % len(lab))
rx = re.compile(r"Predicted structure for (\S+) with length (\d+), pLDDT ([\d.]+), pTM ([\d.]+)")
l31 = sorted(glob.glob(W + "/logs/step31b_esmfold.*.log"), key=os.path.getmtime)[-1]
conf = {}
for lg, want in ((l31, "no"), (W + "/logs/step29b_esmfold.29932169.log", "yes")):
    for m in rx.finditer(open(lg).read()):
        r = m.group(1)
        if r in lab and lab[r]["step29_pdb"] == want: conf[r] = (float(m.group(3)), float(m.group(4)), os.path.basename(lg))
print("ESMFold logs: %s, step29b_esmfold.29932169.log | pLDDT/pTM found for %d of 84" % (os.path.basename(l31), len(conf)))
if len(conf) != 84: raise SystemExit("ABORT: missing pLDDT/pTM for %s" % sorted(set(lab) - set(conf))[:5])
def qid(x):  # as step29e
    m = [q for q in lab if x == q or x.startswith(q + ".") or x.startswith(q + "_")]
    if len(m) != 1: raise SystemExit("ABORT: cannot map query name %s" % x)
    return m[0]
best, nh = {}, {}
for db in ("PDB", "CATH50"):
    for l in open(R + "step31_foldseek_%s.m8" % db):
        f = l.rstrip("\n").split("\t"); q = qid(f[0]); nh[(q, db)] = nh.get((q, db), 0) + 1
        if (q, db) not in best or float(f[8]) > float(best[(q, db)][8]): best[(q, db)] = f
cnt = {}
with open(OUT, "w") as o:
    o.write("rep\tlength\tn_animals\teuk_flag\tdefence_neighbours\tstep29_pdb\tplddt\tptm\tplddt_ca\thelix_frac\tstrand_frac\trg_ratio\tclass")
    for db in ("PDB", "CATH50"):
        o.write("".join("\t%s_%s" % (db, c) for c in ("n_hits", "target", "theader", "prob", "evalue", "qtm", "ttm", "qcov")))
    o.write("\n")
    for r in lab:
        p = Q + "/" + r + ".pdb"
        ca = [l for l in open(p) if l.startswith("ATOM") and l[12:16].strip() == "CA"]
        X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in ca]); bca = np.mean([float(l[60:66]) for l in ca])
        n = len(X); d2 = np.linalg.norm(X[2:] - X[:-2], axis=1); d3 = np.linalg.norm(X[3:] - X[:-3], axis=1)
        hel = float(np.mean((d3 >= 4.5) & (d3 <= 6.0))); strd = float(np.mean((d2[:-1] >= 6.2) & (d3 >= 9.0)))
        rgr = float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean())) / (2.2 * n ** 0.38)
        pl, pt, _ = conf[r]
        cls = ("helix_only" if hel >= 0.8 and strd < 0.05 else "folded") if (pl >= 70 and pt >= 0.5) else "low_confidence"
        cnt[cls] = cnt.get(cls, 0) + 1
        L = lab[r]
        if n != int(L["length"]): raise SystemExit("ABORT: %s has %d CA, length %s" % (r, n, L["length"]))
        row = [r, L["length"], L["n_animals"], L["euk_flag"], L["defence_neighbours(member_frac_pm3)"], L["step29_pdb"],
               "%.1f" % pl, "%.3f" % pt, "%.1f" % bca, "%.2f" % hel, "%.2f" % strd, "%.2f" % rgr, cls]
        for db in ("PDB", "CATH50"):
            b = best.get((r, db))
            row += [str(nh.get((r, db), 0))] + ([b[1], b[2], b[9], b[7], b[11], b[12], b[14]] if b else ["none"] + ["NA"] * 6)
        o.write("\t".join(row) + "\n")
print("classes:", cnt)
print("wrote", OUT)
