# Step 34b2: catalog statistics on strict_final, strict beside it; and the NMPFams2 query filter applied to the
#   strict_final EU reps. Step 34b files are not touched.
# Part 1. Sets: strict = origin Tiara-prokaryotic (results/step34a_rep_origin_v2_final.tsv, column origin);
#   strict_final = column strict_final = 1. Per class (K, KWP, GU, EU, all): reps and % of the set; mean CDS length
#   (bp, fourway_classes_n88.tsv); % complete ORFs (completeness "00"); mean prevalence of 88 samples at >= 1 read and
#   at >= 5 reads; % present in all three cohorts at both criteria; share of mapped reads raw and length-normalised
#   (reads per kb of CDS), pooled and mean of per-sample shares; denominator all reads on the 2,067,011 supported
#   reps (as step34b). Strict at >= 1 read is checked against results/step34b_catalog_stats.tsv.
# Part 2. NMPFams2 query set results/unk_filtered_aplakidou.faa was written by scripts/step26a_unk_linclust30.sh from
#   results/supported_U_reps_n88.faa (585,347 = all GU + EU supported reps) with
#     F1 keep if Prodigal partial=00, or partial=10 (left edge only) with start coordinate > 10 nt
#     F2 keep if length >= 35 aa, terminal "*" not counted.
#   Applied here with the fourway completeness code (Prodigal partial flag) and the rep sequence length
#   (catalog/db/LsPS_AA_95_rep.fasta, "*" stripped); the F1 left-edge exception kept no rep in step26a (all 112,903
#   left-edge partials start within 10 nt), so F1 = partial 00. Prediction is checked against actual membership.
# Outputs (refuse to overwrite): results/step34b2_catalog_stats_final.tsv, results/step34b2_nmpfams2_filter.tsv
import os, gzip, collections, statistics
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
O1, O2 = R + "step34b2_catalog_stats_final.tsv", R + "step34b2_nmpfams2_filter.tsv"
for p in (O1, O2):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
CL = ["K", "KWP", "GU", "EU"]
mem = {}
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); io, ifn = h.index("origin"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        s = f[io] == "Tiara-prokaryotic"; fin = f[ifn] == "1"
        if fin and not s: raise SystemExit("ABORT: strict_final rep outside strict")
        mem[f[0]] = (s, fin)
info = {}
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); info[f[0]] = (f[1], int(f[2]), f[3])
if set(mem) != set(info) or len(info) != 2067011: raise SystemExit("ABORT: tables differ")
coh = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: f = l.split("\t"); coh[f[0]] = f[1]
SETS = ("strict_final", "strict")
G = [(s, c) for s in SETS for c in CL]
st = {g: [0] * 7 for g in G}  # n, len, complete, prev1, prev5, coh3_1, coh3_5
with gzip.open(R + "count_matrix_primary_n88.tsv.gz", "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")[1:]; ns = len(hdr); ci = [coh[s] for s in hdr]
    raw = {g: [0] * ns for g in G}; rpk = {g: [0.0] * ns for g in G}; TR = [0] * ns; TQ = [0.0] * ns
    for l in fh:
        p = l.rstrip("\n").split("\t"); r = p[0]; c, ln, cp = info[r]; kb = ln / 1000.0
        x = [int(v) for v in p[1:]]
        for j in range(ns):
            if x[j]: TR[j] += x[j]; TQ[j] += x[j] / kb
        s, fin = mem[r]
        if not s: continue
        p1 = sum(1 for v in x if v >= 1); p5 = sum(1 for v in x if v >= 5)
        c1 = len({ci[j] for j in range(ns) if x[j] >= 1}) == 3; c5 = len({ci[j] for j in range(ns) if x[j] >= 5}) == 3
        for sn in (("strict_final", "strict") if fin else ("strict",)):
            g = (sn, c); t = st[g]
            t[0] += 1; t[1] += ln; t[2] += cp == "00"; t[3] += p1; t[4] += p5; t[5] += c1; t[6] += c5
            a = raw[g]; b = rpk[g]
            for j in range(ns):
                if x[j]: a[j] += x[j]; b[j] += x[j] / kb
TRs, TQs = sum(TR), sum(TQ)
def row(sn, c):
    gs = [(sn, cc) for cc in (CL if c == "all" else [c])]
    m = sum(st[g][0] for g in gs); tot = sum(st[(sn, cc)][0] for cc in CL)
    S = lambda i: sum(st[g][i] for g in gs)
    rs = [sum(raw[g][j] for g in gs) for j in range(ns)]; qs = [sum(rpk[g][j] for g in gs) for j in range(ns)]
    return [sn, c, m, "%.2f" % (100.0 * m / tot), "%.0f" % (S(1) / m), "%.2f" % (100.0 * S(2) / m), "%.2f" % (S(3) / m), "%.2f" % (S(4) / m),
            "%.2f" % (100.0 * S(5) / m), "%.2f" % (100.0 * S(6) / m), "%.2f" % (100.0 * sum(rs) / TRs), "%.2f" % (100.0 * sum(qs) / TQs),
            "%.2f" % statistics.mean(100.0 * rs[j] / TR[j] for j in range(ns)), "%.2f" % statistics.mean(100.0 * qs[j] / TQ[j] for j in range(ns))]
H = ["set", "class", "n_reps", "pct_of_set", "mean_cds_bp", "pct_complete_orf", "mean_prev_ge1read", "mean_prev_ge5reads",
     "pct_3cohorts_ge1read", "pct_3cohorts_ge5reads", "pct_reads_raw_pooled", "pct_reads_lennorm_pooled",
     "pct_reads_raw_mean_per_sample", "pct_reads_lennorm_mean_per_sample"]
rows = [row(sn, c) for c in CL + ["all"] for sn in SETS]
old = {}
for i, l in enumerate(open(R + "step34b_catalog_stats.tsv")):
    if i: f = l.rstrip("\n").split("\t"); old[(f[0], f[1])] = f
bad = 0
for r in rows:
    if r[0] == "strict":
        o = old[("strict", r[1])]
        if [str(r[2]), r[4], r[5], r[6], r[8], r[10], r[11]] != [o[2], o[4], o[5], o[6], o[7], o[8], o[9]]: bad += 1; print("  differs from step34b:", r[1], o)
print("strict rows reproduced from step34b_catalog_stats.tsv (>= 1 read):", "yes" if not bad else "NO (%d)" % bad, flush=True)
with open(O1, "w") as o:
    o.write("\t".join(H) + "\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("\n%-12s %-4s %9s %6s %6s %6s %13s %13s %13s %13s %13s" % ("set", "cls", "reps", "%set", "bp", "cmpl%", "prev >=1|>=5", "3coh% >=1|>=5", "raw% pool|mean", "lnorm% pool|mean", ""))
for r in rows:
    print("%-12s %-4s %9d %6s %6s %6s %6s|%-6s %6s|%-6s %6s|%-6s %6s|%-6s" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[12], r[11], r[13]))
# Part 2
EU = sorted(r for r, (s, fin) in mem.items() if fin and info[r][0] == "EU")
U = {l[1:].split()[0] for l in open(R + "supported_U_reps_n88.faa") if l[0] == ">"}
Q = {l[1:].split()[0] for l in open(R + "unk_filtered_aplakidou.faa") if l[0] == ">"}
want = set(EU); L = {}; hd = None; n = 0
with open(W + "/catalog/db/LsPS_AA_95_rep.fasta") as fh:
    for l in fh:
        if l[0] == ">":
            if hd in want: L[hd] = n
            hd = l[1:].split()[0]; n = 0
        elif hd in want: n += len(l.strip().replace("*", ""))
    if hd in want: L[hd] = n
print("\nPart 2: strict_final EU reps %d | in supported_U_reps_n88.faa %d | in unk_filtered_aplakidou.faa %d" % (len(EU), len(want & U), len(want & Q)))
cnt = collections.Counter(); mis = 0; codes = collections.Counter()
with open(O2, "w") as o:
    o.write("rep\tpartial_code\tlength_aa\tF1_complete\tF2_ge35aa\tpredicted_in_query_set\tin_query_set\n")
    for r in EU:
        cp = info[r][2]; f1 = cp == "00"; f2 = L[r] >= 35; pr = f1 and f2 and r in U; a = r in Q
        mis += pr != a; codes[cp] += 1
        cnt["fail F1"] += not f1; cnt["fail F2"] += not f2; cnt["fail both"] += (not f1 and not f2)
        cnt["fail F1 only"] += (not f1 and f2); cnt["fail F2 only"] += (f1 and not f2); cnt["pass both"] += f1 and f2
        o.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (r, cp, L[r], f1, f2, pr, a))
print("  partial codes:", dict(codes))
for k in ("pass both", "fail F1", "fail F2", "fail both", "fail F1 only", "fail F2 only"):
    print("  %-13s %6d (%.2f%%)" % (k, cnt[k], 100.0 * cnt[k] / len(EU)))
print("  predicted vs actual query-set membership mismatches: %d" % mis)
print("wrote", O1, "and", O2)
