import glob, os, time, statistics as st
P = "/bigdata/stajichlab/lshad003/wf_protein_space/results/esmfold_pdb"
now = time.time(); rows = []; skipped = 0
for f in sorted(glob.glob(P + "/*.pdb")):
    if now - os.path.getmtime(f) < 120: skipped += 1; continue
    ca = {}
    for l in open(f):
        if l.startswith("ATOM") and l[12:16].strip() == "CA": ca[(l[21], l[22:27])] = float(l[60:66])
    if ca: rows.append([os.path.basename(f)[:-4], list(ca.values())])
if not rows: raise SystemExit("no finished PDB files read")
scale = 1.0 if max(max(v) for _, v in rows) > 1.5 else 100.0
T = []
for name, v in rows:
    v = [x * scale for x in v]; T.append((name, len(v), st.mean(v), sum(x >= 70 for x in v) / len(v)))
L = sorted(t[1] for t in T); M = [t[2] for t in T]
print("structures read %d (skipped %d recently written) | pLDDT scale in file: %s" % (len(T), skipped, "0-100" if scale == 1 else "0-1"))
print("length aa: min %d, median %d, max %d" % (L[0], L[len(L) // 2], L[-1]))
print("mean pLDDT: median %.1f | >=70 %d (%.1f%%) | 50-70 %d | <50 %d" % (st.median(M), sum(m >= 70 for m in M), 100 * sum(m >= 70 for m in M) / len(M),
      sum(50 <= m < 70 for m in M), sum(m < 50 for m in M)))
T.sort(key=lambda t: -t[2])
print("top 10:"); [print("  %-45s len %4d  pLDDT %.1f  frac>=70 %.2f" % t) for t in T[:10]]
print("bottom 5:"); [print("  %-45s len %4d  pLDDT %.1f  frac>=70 %.2f" % t) for t in T[-5:]]
with open("/bigdata/stajichlab/lshad003/wf_protein_space/results/esmfold_plddt_peek.tsv", "w") as o:
    o.write("structure\tlength\tmean_plddt\tfrac_res_ge70\n"); [o.write("%s\t%d\t%.2f\t%.3f\n" % t) for t in T]
