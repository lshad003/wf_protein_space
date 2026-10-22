# Step 34t2: tank vs sequencing run for the Fig. 3f tank contrasts, with final p-values and saved null distributions.
# Units, presence and restriction rule exactly as scripts/step34t_eggmass_tank.py (41,114 strict_final EU, 41,114 random K
#   reps, random.seed(343); assembly presence or >= 5 reads at c_id95; restricted = present in >= 50% of the tank's
#   froglets and >= 3, absent from every froglet of the other tank of the same egg mass; both tanks summed; D = EU - K).
# Run: results/batch_table_88.tsv column run. Nested = every run holds froglets of one tank only (run determines tank).
# WF24 EM4 (13 vs 19, two runs): (a) unrestricted permutation as Step 34t (reproduction: EU 822, K 822 assembly; 1,064,
#   866 c_id95; p as Step 34t), (b) within-run permutation: tank labels permuted within each run, 1,000 permutations,
#   numpy default_rng(34). Final p = (b).
# WF23 EM3 (3 vs 6): every one of the C(9,3) = 84 tank labelings enumerated (the observed one included); exact
#   p = #(labelings with null >= observed) / 84, smallest attainable 1/84 = 0.0119. Final p = exact.
# Permutation p = (1 + #null >= observed) / (n + 1), one-sided. Null summaries: mean, sd (population), 95th percentile
#   (numpy linear interpolation).
# Outputs (refuse to overwrite): results/step34t/tank_run_crosstab.tsv, results/step34t/tank_run_check.tsv,
#   results/step34t/null_<contrast>_<presence>_<method>.tsv (one row per permutation or labeling: EU, K, D).
import os, gzip, random, itertools, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; O = R + "step34t/"
OX, OC = O + "tank_run_crosstab.tsv", O + "tank_run_check.tsv"
for p in (OX, OC):
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
C = collections.OrderedDict()
for nm, emc in (("WF24_EM4_tank", "WF24_EM4"), ("WF23_EM3_tank", "WF23_EM3")):
    C[nm] = {s.split(".")[0]: (tank[s], [s]) for s in stems if meta[s]["cohort"] != "WF22" and em2[s] == emc}
xt = []; nested = {}
for nm, D in C.items():
    c = collections.Counter((D[a][0], run[D[a][1][0]]) for a in D)
    runs = sorted({r for _, r in c}); tanks = sorted({t for t, _ in c})
    nested[nm] = all(sum(1 for t in tanks if c[(t, r)]) == 1 for r in runs) and len(runs) > 1
    print("%s: runs %s; nested %s" % (nm, runs, nested[nm]))
    for t in tanks:
        print("   tank %s: %s" % (t, ", ".join("%s %d" % (r, c[(t, r)]) for r in runs)))
        for r in runs: xt.append((nm, t, r, c[(t, r)]))
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
si = {s: j for j, s in enumerate(stems)}
S = {(nm, m): np.zeros((len(L), len(si)), dtype=bool) for nm, L in SETS.items() for m in ("assembly", "reads>=5_c_id95")}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t"); x = idx.get(r)
        if x is None: continue
        j = si.get(m.split("__", 1)[0])
        if j is not None: S[(x[0], "assembly")][x[1], j] = True
seen = 0
with gzip.open(R + "step34r/matrix_c_id95_n88.tsv.gz", "rt") as fh:
    cols = next(fh).rstrip("\n").split("\t")[1:]
    cj = [(j, si[s]) for j, s in enumerate(cols) if s in si]
    for l in fh:
        r, rest = l.split("\t", 1); x = idx.get(r)
        if x is None: continue
        seen += 1; v = rest.rstrip("\n").split("\t")
        for j, k in cj:
            if float(v[j]) >= 5: S[(x[0], "reads>=5_c_id95")][x[1], k] = True
if seen != len(idx) or len(cj) != len(si): raise SystemExit("ABORT: matrix coverage")
def restricted_total(M, lab):
    tot = M.sum(1); n_ = 0
    for e in sorted(set(lab)):
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1)
        n_ += int(((n >= 3) & (n >= 0.5 * g.sum()) & (tot - n == 0)).sum())
    return n_
labelings = {}
D = C["WF24_EM4_tank"]; an = sorted(D); lab0 = [D[a][0] for a in an]
rng = np.random.default_rng(34); labelings[("WF24_EM4_tank", "unrestricted")] = (lab0, [list(rng.permutation(lab0)) for _ in range(NPERM)])
ra = np.array([run[D[a][1][0]] for a in an]); rng = np.random.default_rng(34); P = []
for _ in range(NPERM):
    lab = np.array(lab0, dtype=object)
    for r in sorted(set(ra)):
        ii = np.where(ra == r)[0]; lab[ii] = rng.permutation(lab[ii])
    P.append(list(lab))
if any(collections.Counter((p[i], ra[i]) for i in range(len(an))) != collections.Counter((lab0[i], ra[i]) for i in range(len(an))) for p in P):
    raise SystemExit("ABORT: within-run permutation changed tank sizes within a run")
labelings[("WF24_EM4_tank", "within_run")] = (lab0, P)
D = C["WF23_EM3_tank"]; an3 = sorted(D); lab3 = [D[a][0] for a in an3]
small = min(set(lab3), key=lab3.count); big = max(set(lab3), key=lab3.count); ns = lab3.count(small)
E = [[small if i in cmb else big for i in range(len(an3))] for cmb in itertools.combinations(range(len(an3)), ns)]
if len(E) != 84 or lab3 not in E: raise SystemExit("ABORT: WF23 enumeration %d" % len(E))
labelings[("WF23_EM3_tank", "exact")] = (lab3, E)
rows = []
for (nm, meth), (l0, L) in labelings.items():
    D = C[nm]; anm = sorted(D)
    for mode in ("assembly", "reads>=5_c_id95"):
        AE = np.stack([S[("EU", mode)][:, [si[s] for s in D[a][1]]].any(1) for a in anm], 1)
        AK = np.stack([S[("K", mode)][:, [si[s] for s in D[a][1]]].any(1) for a in anm], 1)
        oE, oK = restricted_total(AE, l0), restricted_total(AK, l0)
        N = np.array([(restricted_total(AE, p), restricted_total(AK, p)) for p in L]); nd = N[:, 0] - N[:, 1]
        with open(O + "null_%s_%s_%s.tsv" % (nm, mode.replace(">=", "ge"), meth), "w") as o:
            o.write("labeling\t%s\tEU\tK\tD_EU_minus_K\n" % ("tank_" + small + "_froglets" if meth == "exact" else "permutation"))
            for i, (p, x) in enumerate(zip(L, N)):
                key = ",".join(a for a, t in zip(anm, p) if t == small) if meth == "exact" else str(i + 1)
                o.write("%d\t%s\t%d\t%d\t%d\n" % (i + 1, key, x[0], x[1], x[0] - x[1]))
        for st, ob, n in (("EU", oE, N[:, 0]), ("K", oK, N[:, 1]), ("D_EU_minus_K", oE - oK, nd)):
            p = (int((n >= ob).sum())) / len(L) if meth == "exact" else (1 + int((n >= ob).sum())) / (len(L) + 1)
            rows.append((nm, mode, st, meth, ob, len(L), len(set(tuple(x) for x in L)), "%.2f" % n.mean(), "%.2f" % n.std(), "%.2f" % np.percentile(n, 95), "%.4f" % p))
            print("%-14s %-16s %-13s %-12s obs %5d null mean %8.2f sd %7.2f p95 %8.2f p %.4f" % (nm, mode, st, meth, ob, n.mean(), n.std(), np.percentile(n, 95), p), flush=True)
chk = {(r[0], r[1], r[2], r[3]): r for r in rows}
bad = []
for key, exp in ((("WF24_EM4_tank", "assembly", "EU", "unrestricted"), (822, "0.0010")), (("WF24_EM4_tank", "assembly", "K", "unrestricted"), (822, "0.0010")),
                 (("WF24_EM4_tank", "assembly", "D_EU_minus_K", "unrestricted"), (0, "0.6194")),
                 (("WF24_EM4_tank", "reads>=5_c_id95", "EU", "unrestricted"), (1064, "0.0010")), (("WF24_EM4_tank", "reads>=5_c_id95", "K", "unrestricted"), (866, "0.0010")),
                 (("WF24_EM4_tank", "reads>=5_c_id95", "D_EU_minus_K", "unrestricted"), (198, "0.0010"))):
    ok = chk[key][4] == exp[0] and chk[key][10] == exp[1]
    print("CHECK %-60s got %s p %s expected %s p %s %s" % (" ".join(key), chk[key][4], chk[key][10], exp[0], exp[1], "PASS" if ok else "FAIL"))
    if not ok: bad.append(" ".join(key))
with open(OX, "w") as o:
    o.write("contrast\ttank\trun\tn_froglets\n"); o.writelines("%s\t%s\t%s\t%d\n" % r for r in xt)
    for nm in C: o.write("%s\tnested\t-\t%d\n" % (nm, int(nested[nm])))
with open(OC, "w") as o:
    o.write("contrast\tpresence\tset\tmethod\tobserved_restricted\tn_null\tdistinct_labelings\tnull_mean\tnull_sd\tnull_p95\tp_one_sided\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("checks failed: %d %s" % (len(bad), bad)); print("wrote", OX, OC)
raise SystemExit(1 if bad else 0)
