# Step 34s reproduction check: every output written from the old matrix (results/step34s/old/, scripts/step34s_run.sh)
#   against the original output in results/. Byte comparison (gzip files decompressed first); step34k_summary.tsv is
#   compared on its Part D rows only (the reproduction runs Part D only). Exits 1 on any difference.
# Output: results/step34s/repro_old_check.tsv; refuses to overwrite.
import os, gzip, hashlib, sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; N = R + "step34s/old/"
OUT = R + "step34s/repro_old_check.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
def md5(p):
    h = hashlib.md5(); f = gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")
    for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
F = ["step34b2_catalog_stats_final.tsv", "step34b2_nmpfams2_filter.tsv", "prevalence_primary_n88.tsv.gz",
     "fig_prep_read_retention.tsv", "step34e_counts_strict_final.tsv.gz", "step34e_counts_strict_final_EU.tsv.gz",
     "step34e_strict_final_permanova.tsv", "step34e_strict_final_limma.tsv", "step34e_strict_final_EU_permanova.tsv",
     "step34e_strict_final_EU_limma.tsv", "step34e_summary.tsv", "step34l_m23_permanova.tsv", "fig_pcoa_wf22.tsv"]
F += ["step34e_limma/" + f for f in sorted(os.listdir(R + "step34e_limma"))]
rows = []; bad = 0
for f in F:
    a, b = R + f, N + f
    if not os.path.exists(b): rows.append((f, md5(a), "MISSING", "0")); bad += 1; continue
    ma, mb = md5(a), md5(b); ok = ma == mb; bad += not ok; rows.append((f, ma, mb, str(int(ok))))
extra = sorted(set(os.listdir(N + "step34e_limma")) - set(os.listdir(R + "step34e_limma")))
for f in extra: rows.append(("step34e_limma/" + f, "ABSENT_IN_OLD", md5(N + "step34e_limma/" + f), "0")); bad += 1
dold = [l for l in open(R + "step34k_summary.tsv") if l.startswith("D\t")]
dnew = [l for l in open(N + "step34k_summary.tsv") if l.startswith("D\t")]
# source_file column differs by design (path of the new retention table); compare section, metric, value
key = lambda ls: [l.split("\t")[:3] for l in ls]
ok = key(dold) == key(dnew) and len(dold) > 0; bad += not ok
rows.append(("step34k_summary.tsv (Part D rows: section, metric, value)", "%d rows" % len(dold), "%d rows" % len(dnew), str(int(ok))))
with open(OUT, "w") as o:
    o.write("file\tmd5_original\tmd5_reproduced\tidentical\n")
    for r in rows: o.write("\t".join(r) + "\n")
for r in rows: print("%s  %s" % ("OK  " if r[3] == "1" else "DIFF", r[0]))
print("files compared %d, differing %d; wrote %s" % (len(rows), bad, OUT))
sys.exit(1 if bad else 0)
