# Step 34t3: tank audit for the Fig. 3f WF24 EM4 contrast (no figure changes).
# 1. Depth and assembly size of the 32 WF24 EM4 froglets (UHM586.41010, UHM590.41012 excluded as Step 34t): post-fastp
#    read pairs = fastp after_filtering total_reads / 2 (results/step34r/fastp/<stem>.json, the Step 34r trimming);
#    assembly size = scaf_bp of the MEGAHIT assembly ($FEC/<stem>/stats_scaffolds.tsv). Median and range by tank within
#    each run (results/batch_table_88.tsv); two-sided Wilcoxon rank-sum (Mann-Whitney U, scipy 1.13.1, exact when no
#    ties) of tank 4.1 vs 4.2, pooled across runs.
# 2. Observed tank-restricted counts split by direction (restricted to 4.1 vs to 4.2), assembly and >= 5 reads at c_id95,
#    rule and units exactly as scripts/step34t_eggmass_tank.py; restricted ID lists saved; EU and K lists checked disjoint.
# 3. Only if the read-pair Wilcoxon p < 0.05: each EM4 sample's c_id95 counts binomially thinned to the smallest EM4
#    c_id95 total (results/step34r/matrix_colsums_n88.tsv, all genes; p_keep = min total / sample total; numpy
#    default_rng(343)), tank restriction recounted at >= 5 reads, within-run test (1,000 permutations, default_rng(34),
#    as scripts/step34t2_tank_run_check.py) for EU, K and EU - K.
# 4. Egg mass x sequencing run for the WF22 months 2-3 samples (metadata/wf22_design.tsv).
# Outputs (refuse to overwrite) in results/step34t/: tank_audit_depth.tsv, tank_audit_depth_summary.tsv,
#   tank_audit_direction.tsv, wf24_tank_restricted_ids.tsv, tank_audit_thinned.tsv (only if item 3 runs),
#   wf22_eggmass_run_crosstab.tsv.
import os, gzip, json, random, collections
import numpy as np
from scipy.stats import mannwhitneyu
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; O = R + "step34t/"
FEC = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results"
OUT = {k: O + v for k, v in (("depth", "tank_audit_depth.tsv"), ("dsum", "tank_audit_depth_summary.tsv"), ("dir", "tank_audit_direction.tsv"),
       ("ids", "wf24_tank_restricted_ids.tsv"), ("thin", "tank_audit_thinned.tsv"), ("x22", "wf22_eggmass_run_crosstab.tsv"))}
for p in OUT.values():
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
NPERM = 1000
meta = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); meta[f["stem"]] = f
excl = {l.split("\t")[0] for l in list(open(W + "/metadata/wf24_excluded.tsv"))[1:]}
if excl != {"UHM586.41010", "UHM590.41012"}: raise SystemExit("ABORT: exclusions %s" % excl)
tank, em2 = {}, {}
with open(W + "/metadata/wf23_wf24_sequenced_metadata.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); tank[f["stem"]] = f["tank"]; em2[f["stem"]] = f["wf_year"] + "_EM" + f["egg_mass"]
run = {}
with open(R + "batch_table_88.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); run[f["stem"]] = f["run"]
stems = sorted(s for s in meta if s not in excl)
em4 = sorted((s for s in stems if meta[s]["cohort"] != "WF22" and em2[s] == "WF24_EM4"), key=lambda s: s.split(".")[0])
if len(em4) != 32 or len({s.split(".")[0] for s in em4}) != 32: raise SystemExit("ABORT: EM4 froglets %d" % len(em4))
# 1
dep = {}
for s in em4:
    rp = json.load(open(R + "step34r/fastp/%s.json" % s))["summary"]["after_filtering"]["total_reads"] / 2
    st = list(open(FEC + "/%s/stats_scaffolds.tsv" % s)); hh = st[0].rstrip("\n").split("\t"); v = st[1].rstrip("\n").split("\t")
    dep[s] = (rp, int(v[hh.index("scaf_bp")]))
with open(OUT["depth"], "w") as o:
    o.write("stem\ttank\trun\tpost_fastp_read_pairs\tassembly_bp\n")
    for s in em4: o.write("%s\t%s\t%s\t%d\t%d\n" % (s, tank[s], run[s], dep[s][0], dep[s][1]))
summ = []; wil = {}
for mi, meas in enumerate(("post_fastp_read_pairs", "assembly_bp")):
    print("\n1. %s" % meas)
    for r in sorted({run[s] for s in em4}):
        for t in ("4.1", "4.2"):
            v = [dep[s][mi] for s in em4 if run[s] == r and tank[s] == t]
            summ.append((meas, r, t, len(v), "%.0f" % np.median(v), "%.0f" % min(v), "%.0f" % max(v), "", ""))
            print("   %-22s tank %s n %2d median %14s range %s - %s" % (r, t, len(v), format(int(np.median(v)), ","), format(int(min(v)), ","), format(int(max(v)), ",")))
    a = [dep[s][mi] for s in em4 if tank[s] == "4.1"]; b = [dep[s][mi] for s in em4 if tank[s] == "4.2"]
    U, p = mannwhitneyu(a, b, alternative="two-sided"); wil[meas] = p
    summ.append((meas, "pooled", "4.1 vs 4.2", "%d vs %d" % (len(a), len(b)), "%.0f vs %.0f" % (np.median(a), np.median(b)),
                 "%.0f-%.0f" % (min(a), max(a)), "%.0f-%.0f" % (min(b), max(b)), "%.1f" % U, "%.4f" % p))
    print("   pooled: 4.1 median %s (n %d) vs 4.2 median %s (n %d); Wilcoxon rank-sum U %.1f, two-sided p %.4f"
          % (format(int(np.median(a)), ","), len(a), format(int(np.median(b)), ","), len(b), U, p))
with open(OUT["dsum"], "w") as o:
    o.write("measure\trun\ttank\tn\tmedian\tmin_or_range41\tmax_or_range42\twilcoxon_U\tp_two_sided\n")
    for r in summ: o.write("\t".join(map(str, r)) + "\n")
# units
EU, K = [], []
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, ifn = h.index("class4"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ifn] == "1":
            if f[ic] == "EU": EU.append(f[0])
            elif f[ic] == "K": K.append(f[0])
if len(EU) != 41114: raise SystemExit("ABORT: expected 41,114 EU")
random.seed(343); Kb = sorted(random.sample(sorted(K), len(EU)))
SETS = {"EU": sorted(EU), "K": Kb}
idx = {r: (nm, i) for nm, L in SETS.items() for i, r in enumerate(L)}
ej = {s: j for j, s in enumerate(em4)}
ASM = {nm: np.zeros((len(L), 32), dtype=bool) for nm, L in SETS.items()}
CNT = {nm: np.zeros((len(L), 32), dtype=np.int64) for nm, L in SETS.items()}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t"); x = idx.get(r)
        if x is None: continue
        j = ej.get(m.split("__", 1)[0])
        if j is not None: ASM[x[0]][x[1], j] = True
seen = 0
with gzip.open(R + "step34r/matrix_c_id95_n88.tsv.gz", "rt") as fh:
    cols = next(fh).rstrip("\n").split("\t")[1:]
    cj = [(j, ej[s]) for j, s in enumerate(cols) if s in ej]
    for l in fh:
        r, rest = l.split("\t", 1); x = idx.get(r)
        if x is None: continue
        seen += 1; v = rest.rstrip("\n").split("\t")
        for j, k in cj: CNT[x[0]][x[1], k] = int(float(v[j]))
if seen != len(idx) or len(cj) != 32: raise SystemExit("ABORT: matrix coverage")
lab0 = [tank[s] for s in em4]; ra = np.array([run[s] for s in em4])
def restricted(M, lab):
    tot = M.sum(1); out = {}
    for e in sorted(set(lab)):
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1); out[e] = (n >= 3) & (n >= 0.5 * g.sum()) & (tot - n == 0)
    return out
# 2
drows = []; ids = []; exp = {("assembly", "EU"): 822, ("assembly", "K"): 822, ("reads>=5_c_id95", "EU"): 1064, ("reads>=5_c_id95", "K"): 866}
print("\n2. direction of the observed WF24 EM4 tank-restricted counts")
for mode in ("assembly", "reads>=5_c_id95"):
    sets_ = {}
    for st in ("EU", "K"):
        M = ASM[st] if mode == "assembly" else CNT[st] >= 5
        rs = restricted(M, lab0); tot = int(rs["4.1"].sum() + rs["4.2"].sum())
        if tot != exp[(mode, st)]: raise SystemExit("ABORT: %s %s total %d != %d" % (mode, st, tot, exp[(mode, st)]))
        if (rs["4.1"] & rs["4.2"]).any(): raise SystemExit("ABORT: rep restricted to both tanks")
        drows.append((mode, st, int(rs["4.1"].sum()), int(rs["4.2"].sum()), tot))
        print("   %-16s %-2s restricted to 4.1 (13 froglets) %4d | to 4.2 (19 froglets) %4d | total %4d (reproduced)" % (mode, st, int(rs["4.1"].sum()), int(rs["4.2"].sum()), tot))
        sets_[st] = set()
        for t in ("4.1", "4.2"):
            for i in np.where(rs[t])[0]:
                ids.append((SETS[st][i], st, mode, t)); sets_[st].add(SETS[st][i])
    ov = sets_["EU"] & sets_["K"]
    print("   %-16s EU and K restricted ID lists disjoint: %s (shared %d)" % (mode, "yes" if not ov else "NO", len(ov)))
    drows.append((mode, "EU_K_shared_ids", len(ov), "", ""))
with open(OUT["dir"], "w") as o:
    o.write("presence\tset\trestricted_to_4.1\trestricted_to_4.2\ttotal\n")
    for r in drows: o.write("\t".join(map(str, r)) + "\n")
with open(OUT["ids"], "w") as o:
    o.write("rep\tset\tpresence\trestricted_to_tank\n"); o.writelines("%s\t%s\t%s\t%s\n" % r for r in ids)
# 3
print("\n3. read-pair Wilcoxon p %.4f" % wil["post_fastp_read_pairs"])
if wil["post_fastp_read_pairs"] < 0.05:
    tot = {}
    with open(R + "step34r/matrix_colsums_n88.tsv") as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh:
            f = dict(zip(h, l.rstrip("\n").split("\t"))); tot[f["sample"]] = int(f["c_id95"])
    mn = min(tot[s] for s in em4); pk = np.array([mn / tot[s] for s in em4])
    print("   thinning to the smallest EM4 c_id95 total %s (%s); keep fraction %.3f-%.3f" % (format(mn, ","), min(em4, key=lambda s: tot[s]), pk.min(), pk.max()))
    rng = np.random.default_rng(343)
    TH = {st: rng.binomial(CNT[st], pk[None, :]) >= 5 for st in ("EU", "K")}
    rng = np.random.default_rng(34); P = []
    for _ in range(NPERM):
        lab = np.array(lab0, dtype=object)
        for r in sorted(set(ra)):
            ii = np.where(ra == r)[0]; lab[ii] = rng.permutation(lab[ii])
        P.append(list(lab))
    sc = lambda lab: tuple(sum(int(v.sum()) for v in restricted(TH[st], lab).values()) for st in ("EU", "K"))
    oE, oK = sc(lab0); N = np.array([sc(p) for p in P]); trows = []
    rs = {st: restricted(TH[st], lab0) for st in ("EU", "K")}
    for st, ob, n in (("EU", oE, N[:, 0]), ("K", oK, N[:, 1]), ("D_EU_minus_K", oE - oK, N[:, 0] - N[:, 1])):
        p = (1 + int((n >= ob).sum())) / (NPERM + 1)
        d41 = int(rs[st]["4.1"].sum()) if st != "D_EU_minus_K" else int(rs["EU"]["4.1"].sum() - rs["K"]["4.1"].sum())
        d42 = int(rs[st]["4.2"].sum()) if st != "D_EU_minus_K" else int(rs["EU"]["4.2"].sum() - rs["K"]["4.2"].sum())
        trows.append((st, ob, d41, d42, "%.2f" % n.mean(), "%.2f" % n.std(), "%.2f" % np.percentile(n, 95), "%.4f" % p))
        print("   thinned >= 5 reads %-13s obs %5d (4.1 %d, 4.2 %d) null mean %.2f sd %.2f p95 %.2f within-run p %.4f" % (st, ob, d41, d42, n.mean(), n.std(), np.percentile(n, 95), p))
    with open(OUT["thin"], "w") as o:
        o.write("set\tobserved_restricted\tto_4.1\tto_4.2\tnull_mean\tnull_sd\tnull_p95\tp_one_sided_within_run\n")
        for r in trows: o.write("\t".join(map(str, r)) + "\n")
else:
    print("   not < 0.05: thinning not run")
# 4
d22 = {}
for l in list(open(W + "/metadata/wf22_design.tsv"))[1:]:
    s, a, t, e, m = l.rstrip("\n").split("\t"); d22[s] = (a, int(m), "EM" + e)
w = [s for s in stems if meta[s]["cohort"] == "WF22" and d22[s][1] in (2, 3)]
c = collections.Counter((d22[s][2], d22[s][1], run[s]) for s in w); runs = sorted({run[s] for s in w})
print("\n4. WF22 months 2-3: egg mass x sequencing run (samples; %d samples, %d froglets)" % (len(w), len({d22[s][0] for s in w})))
print("   %-12s %s" % ("", "  ".join("%22s" % r for r in runs)))
with open(OUT["x22"], "w") as o:
    o.write("egg_mass\tmonth\t" + "\t".join(runs) + "\n")
    for e in ("EM1", "EM2", "EM3"):
        for m in (2, 3):
            o.write("%s\t%d\t%s\n" % (e, m, "\t".join(str(c[(e, m, r)]) for r in runs)))
            print("   %s month %d %s" % (e, m, "  ".join("%22d" % c[(e, m, r)] for r in runs)))
print("\nwrote", ", ".join(p for p in OUT.values() if os.path.exists(p)))
