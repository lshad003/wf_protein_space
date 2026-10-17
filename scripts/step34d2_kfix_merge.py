# Step 34d2 merge: K hit rates on the 4,509 distinct reps of results/step34d_known_lenmatched_fixed.faa.
# Hits: Step 34d rows (results/step34d_prostt5/known_vs_{PDB,CATH50}.tsv) for the 4,284 retained reps plus the
#   225 replacement reps (results/step34d2_prostt5/krepl_vs_{PDB,CATH50}.tsv), same settings. Hit = best E <= 1e-3 in
#   either database; bins as step34d_summary.py. Query and decoy rows copied from results/step34d_by_length.tsv.
# Output: results/step34d2_by_length_fixed.tsv (old file kept); refuses to overwrite.
import os
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = R + "step34d2_by_length_fixed.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
BINS = (("<100", 0, 100), ("100-199", 100, 200), ("200-399", 200, 400), (">=400", 400, 10 ** 9), ("all", 0, 10 ** 9))
L = {}; h = None
for l in open(R + "step34d_known_lenmatched_fixed.faa"):
    if l[0] == ">": h = l[1:].strip()
    else: L[h] = len(l.strip())
new = {l[1:].strip() for l in open(R + "step34d2_k_replacements.faa") if l[0] == ">"}
if len(L) != 4509 or len(new) != 225 or not new <= set(L): raise SystemExit("ABORT: fixed set")
b = {}
for f, ids in ((R + "step34d_prostt5/known_vs_%s.tsv", set(L) - new), (R + "step34d2_prostt5/krepl_vs_%s.tsv", new)):
    for db in ("PDB", "CATH50"):
        for l in open(f % db):
            x = l.split("\t")
            if x[0] not in ids: continue
            e = float(x[2]); b[x[0]] = min(e, b.get(x[0], 9e9))
H = {k for k, e in b.items() if e <= 1e-3}
print("replacements: with any row %d, hit %d of 225" % (sum(1 for k in new if k in b), len(H & new)))
rows = [l.rstrip("\n").split("\t") for l in open(R + "step34d_by_length.tsv")]
out = [rows[0]]
for r in rows[1:]:
    if r[0] != "K_lenmatched": out.append(r); continue
    lo, hi = [(lo, hi) for bn, lo, hi in BINS if bn == r[1]][0]
    ids = [k for k in L if lo <= L[k] < hi]; n = len(ids); k = sum(1 for x in ids if x in H)
    print("K %-8s old %s/%s (%s%%) -> fixed %d/%d (%.2f%%)" % (r[1], r[3], r[2], r[4], k, n, 100.0 * k / n))
    out.append(["K_lenmatched", r[1], str(n), str(k), "%.2f" % (100.0 * k / n)])
with open(OUT, "w") as o:
    for r in out: o.write("\t".join(r) + "\n")
print("wrote", OUT)
