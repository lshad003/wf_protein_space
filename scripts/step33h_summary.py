# Step 33h summary: classes of results/step33h_rep_class.tsv with DIAMOND blastp (--sensitive) hits against the dm6
# proteome (results/step33h_diamond_dm6.tsv, outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore, best hit).
# Hit = any row with evalue <= 1e-3 (also reported at <= 1e-10). Per set (1: EU supported reps on fly95 contigs;
# 2: reps of the fly-majority families whose own contig is fly95) and per class: n, fraction, soft-masked >= 50% of the
# projected span, mobile_genetic_element overlap, DIAMOND hit rate and median identity.
# Output (refuses to overwrite): results/step33h_summary.tsv
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = R + "step33h_summary.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
CL = ["CDS same strand, in frame", "CDS same strand, other frame", "CDS opposite strand", "UTR or non-coding exon", "intron", "intergenic", "not projected"]
rows = [l.rstrip("\n").split("\t") for i, l in enumerate(open(R + "step33h_rep_class.tsv")) if i]
dm = {}
for l in open(R + "step33h_diamond_dm6.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] not in dm or float(f[7]) > float(dm[f[0]][7]): dm[f[0]] = f
fam_total = sum(1 for i, l in enumerate(open(R + "step33g_family_fly.tsv")) if i and l.rstrip("\n").split("\t")[4] == "1")
out = [("set", "class", "n", "pct_of_set", "n_softmasked_ge50", "pct_softmasked_ge50", "n_mge", "n_diamond_e1e-3", "pct_diamond_e1e-3", "n_diamond_e1e-10", "median_pident_e1e-3")]
for si, lab in ((1, "set1_fly_EU_reps"), (2, "set2_flymajority_family_reps")):
    S = [r for r in rows if r[si] == "1"]; n = len(S)
    print("\n%s: n=%d%s" % (lab, n, "" if si == 1 else " (of %d fly-majority families; %d family reps are not on a fly95 contig)" % (fam_total, fam_total - n)))
    print("  %-30s %8s %8s %12s %6s %16s %16s %8s" % ("class", "n", "%", "softmask>=50%", "MGE", "DIAMOND e<=1e-3", "e<=1e-10", "med pid"))
    for k in CL + ["all"]:
        T = S if k == "all" else [r for r in S if r[3] == k]; m = len(T)
        if m == 0: continue
        sm = sum(1 for r in T if float(r[8]) >= 0.5); te = sum(1 for r in T if r[9] == "1")
        h3 = [dm[r[0]] for r in T if r[0] in dm and float(dm[r[0]][6]) <= 1e-3]; h10 = sum(1 for r in T if r[0] in dm and float(dm[r[0]][6]) <= 1e-10)
        pid = sorted(float(x[2]) for x in h3); mp = "%.1f" % pid[len(pid) // 2] if pid else "NA"
        out.append((lab, k, m, "%.2f" % (100.0 * m / n), sm, "%.2f" % (100.0 * sm / m), te, len(h3), "%.2f" % (100.0 * len(h3) / m), h10, mp))
        print("  %-30s %8d %8.2f %6d (%4.1f%%) %6d %7d (%5.1f%%) %7d (%5.1f%%) %8s" % (k, m, 100.0 * m / n, sm, 100.0 * sm / m, te, len(h3), 100.0 * len(h3) / m, h10, 100.0 * h10 / m, mp))
with open(OUT, "w") as o:
    for r in out: o.write("\t".join(map(str, r)) + "\n")
print("\nwrote", OUT)
