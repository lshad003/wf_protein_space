# Step 34b: catalog statistics on the strict set (primary) and the broad set (sensitivity), with the old whole-catalog
# n88 values beside them.
# Inputs (printed at run time):
#   results/step34a_rep_origin.tsv            origin per supported rep (step34a)
#   results/fourway_classes_n88.tsv           class4, CDS length (bp, from supported_cds_lengths.tsv), completeness
#                                             ("00" = complete ORF), n_samples, n_years (step13i_fourway_n88.py)
#   results/count_matrix_primary_n88.tsv.gz   reads per supported rep x 88 samples (2,067,011 rows)
#   metadata/samples_88_treatment.tsv         cohort per sample
# Sets: whole = all 2,067,011 supported reps; strict = origin Tiara-prokaryotic (contigs >= 1 kb, not fly or frog);
#   broad = strict + Tiara-unknown + <1 kb unclassified (not fly or frog).
# Per set x class (K, KWP, GU, EU, all): reps and % of the set; mean CDS length; % complete ORFs; mean prevalence
#   (samples with >= 1 read, from the count matrix; cross-checked against n_samples); % present in all three cohorts
#   (cross-checked against n_years); share of mapped reads, raw and length-normalised (reads per kb of CDS, as
#   step16w), as a pooled share over all 88 samples and as the mean of per-sample shares. Denominator for every share:
#   all reads on the 2,067,011 supported reps (as step13i/step16w).
# Diet: per sample, share of mapped reads on fly-origin reps (also frog and Tiara-eukaryotic other), raw and
#   length-normalised, summarised by cohort.
# Outputs (refuse to overwrite): results/step34b_catalog_stats.tsv, results/step34b_fly_read_share.tsv
import os, gzip, collections, statistics
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = [R + "step34b_catalog_stats.tsv", R + "step34b_fly_read_share.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
IN = [R + "step34a_rep_origin.tsv", R + "fourway_classes_n88.tsv", R + "count_matrix_primary_n88.tsv.gz", W + "/metadata/samples_88_treatment.tsv"]
print("inputs:\n  " + "\n  ".join(IN), flush=True)
ORI = ["fly", "frog", "Tiara-eukaryotic other", "Tiara-prokaryotic", "Tiara-unknown", "<1 kb unclassified"]
CL = ["K", "KWP", "GU", "EU"]
SETS = (("strict", ["Tiara-prokaryotic"]), ("broad", ["Tiara-prokaryotic", "Tiara-unknown", "<1 kb unclassified"]), ("whole_n88", ORI))
ori = {}
for i, l in enumerate(open(IN[0])):
    if i: f = l.rstrip("\n").split("\t"); ori[f[0]] = f[2]
info = {}
for i, l in enumerate(open(IN[1])):
    if i: f = l.rstrip("\n").split("\t"); info[f[0]] = (f[1], int(f[2]), f[3] == "00", int(f[4]), int(f[5]))
if len(ori) != 2067011 or set(ori) != set(info): raise SystemExit("ABORT: origin and fourway tables do not cover the same 2,067,011 reps")
coh = {}
for i, l in enumerate(open(IN[3])):
    if i: f = l.split("\t"); coh[f[0]] = f[1]
G = [(o, c) for o in ORI for c in CL]
st = {g: [0, 0, 0, 0, 0] for g in G}   # n, sum length, n complete, sum prevalence, n in 3 cohorts
mism = collections.Counter()
with gzip.open(IN[2], "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")[1:]
    if len(hdr) != 88 or set(hdr) != set(coh): raise SystemExit("ABORT: count matrix columns are not the 88 samples")
    ns = len(hdr); ci = [coh[s] for s in hdr]
    raw = {g: [0] * ns for g in G}; rpk = {g: [0.0] * ns for g in G}; n = 0
    for l in fh:
        p = l.rstrip("\n").split("\t"); r = p[0]
        c, ln, comp, nsamp, nyr = info[r]; g = (ori[r], c); kb = ln / 1000.0
        a = raw[g]; b = rpk[g]; prev = 0; cs = set()
        for j in range(ns):
            x = int(p[j + 1])
            if x: a[j] += x; b[j] += x / kb; prev += 1; cs.add(ci[j])
        s = st[g]; s[0] += 1; s[1] += ln; s[2] += comp; s[3] += prev; s[4] += len(cs) == 3
        mism["prevalence"] += prev != nsamp; mism["cohorts"] += len(cs) != nyr
        n += 1
        if n % 500000 == 0: print("rows", n, flush=True)
if n != 2067011: raise SystemExit("ABORT: %d count-matrix rows" % n)
print("count matrix rows %d | mismatches vs fourway table: prevalence %d, cohorts %d" % (n, mism["prevalence"], mism["cohorts"]), flush=True)
TR = [sum(raw[g][j] for g in G) for j in range(ns)]; TQ = [sum(rpk[g][j] for g in G) for j in range(ns)]
TRs, TQs = sum(TR), sum(TQ)
rows = [("set", "class", "n_reps", "pct_of_set", "mean_cds_bp", "pct_complete_orf", "mean_prevalence_of_88", "pct_in_3_cohorts",
         "pct_reads_raw_pooled", "pct_reads_lennorm_pooled", "pct_reads_raw_mean_per_sample", "pct_reads_lennorm_mean_per_sample")]
for sn, oris in SETS:
    tot_n = sum(st[(o, c)][0] for o in oris for c in CL)
    for c in CL + ["all"]:
        gs = [(o, cc) for o in oris for cc in (CL if c == "all" else [c])]
        m = sum(st[g][0] for g in gs)
        if m == 0: continue
        L = sum(st[g][1] for g in gs) / m; cp = 100.0 * sum(st[g][2] for g in gs) / m; pv = sum(st[g][3] for g in gs) / m
        y3 = 100.0 * sum(st[g][4] for g in gs) / m
        rs = [sum(raw[g][j] for g in gs) for j in range(ns)]; qs = [sum(rpk[g][j] for g in gs) for j in range(ns)]
        rows.append((sn, c, m, "%.2f" % (100.0 * m / tot_n), "%.0f" % L, "%.2f" % cp, "%.2f" % pv, "%.2f" % y3,
                     "%.2f" % (100.0 * sum(rs) / TRs), "%.2f" % (100.0 * sum(qs) / TQs),
                     "%.2f" % statistics.mean(100.0 * rs[j] / TR[j] for j in range(ns)), "%.2f" % statistics.mean(100.0 * qs[j] / TQ[j] for j in range(ns))))
with open(OUT[0], "w") as o:
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("\nper set x class (shares: % of reads on all 2,067,011 supported reps; pooled | mean per sample)")
print("%-10s %-4s %10s %7s %8s %9s %9s %8s %15s %15s" % ("set", "cls", "reps", "%set", "mean_bp", "complete%", "mean_prev", "3coh%", "raw% pool|mean", "lennorm% pool|mean"))
for r in rows[1:]:
    print("%-10s %-4s %10d %7s %8s %9s %9s %8s %7s|%-7s %7s|%-7s" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[10], r[9], r[11]))
with open(OUT[1], "w") as o:
    o.write("stem\tcohort\tpct_reads_fly_raw\tpct_reads_fly_lennorm\tpct_reads_frog_raw\tpct_reads_tiara_euk_other_raw\n")
    per = []
    for j, s in enumerate(hdr):
        fr = sum(raw[("fly", c)][j] for c in CL); fq = sum(rpk[("fly", c)][j] for c in CL)
        gr = sum(raw[("frog", c)][j] for c in CL); er = sum(raw[("Tiara-eukaryotic other", c)][j] for c in CL)
        v = (100.0 * fr / TR[j], 100.0 * fq / TQ[j], 100.0 * gr / TR[j], 100.0 * er / TR[j]); per.append((s, coh[s]) + v)
        o.write("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n" % ((s, coh[s]) + v))
print("\nper-sample share of mapped reads on fly-origin reps, by cohort (descriptive; cohort is confounded with run)")
for c in ("WF22", "WF23", "WF24", "all"):
    P = [x for x in per if c == "all" or x[1] == c]
    for i, lab in ((2, "fly raw"), (3, "fly length-normalised"), (4, "frog raw"), (5, "Tiara-euk other raw")):
        v = sorted(x[i] for x in P)
        print("  %-4s n=%2d %-22s median %6.2f%%  mean %6.2f%%  range %6.2f%% - %6.2f%%" % (c, len(P), lab, statistics.median(v), statistics.mean(v), v[0], v[-1]))
for p in OUT: print("wrote", p)
