# Step 34m: does rearing tank alone produce restricted unknown proteins within one clutch? (descriptive, permutation null)
# Contrasts: WF24 egg mass 4, tank 4.1 vs 4.2; WF23 egg mass 3, tank 5.0 vs 6.0 (tank from
#   metadata/wf23_wf24_sequenced_metadata.tsv). WF22 egg-mass restriction repeated with one sample per animal (month 2,
#   and separately month 3; metadata/wf22_design.tsv) so WF22 and WF24 both use one sample per animal; WF22 months 2-3
#   pooled is rerun as a reproduction check against Step 34c3b (EU 948, K 616, D 332, null D -7.1, assembly).
# Exclusions: UHM586.41010 and UHM590.41012 (metadata/wf24_excluded.tsv); UHM586 is already absent from the n88 set.
# Units and presence as scripts/step34c3_eggmass.py: 41,114 strict_final EU reps and 41,114 random strict_final K reps
#   (seed 343); presence in an animal = any assembled member (primary, results/clusters_95_n88.tsv) or >= 5 reads in any
#   of its retained samples (sensitivity, results/count_matrix_primary_n88.tsv.gz).
# Restricted (group g within a contrast): present in >= 50% of g's animals AND >= 3, absent from every other animal of
#   the contrast (the other tank of the same egg mass; the other WF22 egg masses).
# Null as scripts/step34c3b_checks.py: group labels permuted among the contrast's animals (sizes kept), 1,000
#   permutations from numpy default_rng(34) (fresh generator per contrast, animals sorted), EU and K scored on the same
#   permutations. p = (1 + #null >= observed) / 1,001, one-sided, for EU, K and D = EU - K; ratio = observed / null mean.
# Overlap: WF24 / WF23 tank-restricted EU reps vs WF22 egg-mass-restricted EU reps (months 2-3 pooled and per month).
# Outputs: results/step34m_tank_restriction.tsv, results/step34m_summary.tsv; refuses to overwrite.
import os, gzip, random, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT1 = R + "step34m_tank_restriction.tsv"; OUT2 = R + "step34m_summary.tsv"
for p in (OUT1, OUT2):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
NPERM = 1000
meta = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); meta[f["stem"]] = f
if len(meta) != 88: raise SystemExit("ABORT: %d stems" % len(meta))
excl = set()
with open(W + "/metadata/wf24_excluded.tsv") as fh:
    next(fh)
    for l in fh: excl.add(l.split("\t")[0])
if excl != {"UHM586.41010", "UHM590.41012"}: raise SystemExit("ABORT: unexpected exclusions %s" % excl)
tank = {}; trt2 = {}; em2 = {}
with open(W + "/metadata/wf23_wf24_sequenced_metadata.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t")))
        tank[f["stem"]] = f["tank"]; trt2[f["stem"]] = f["treatment"]; em2[f["stem"]] = f["wf_year"] + "_EM" + f["egg_mass"]
d22 = {}
with open(W + "/metadata/wf22_design.tsv") as fh:
    next(fh)
    for l in fh:
        s, a, t, e, m = l.rstrip("\n").split("\t"); d22[s] = (a, int(m), t, "WF22_EM" + e)
summ = []; checks = []
def check(name, ok, note=""):
    checks.append((name, "PASS" if ok else "FAIL", note)); print("CHECK %-45s %s %s" % (name, "PASS" if ok else "FAIL", note))
stems = [s for s in meta if s not in excl]
check("UHM586 absent from n88", "UHM586.41010" not in meta)
check("UHM590 present in n88 and removed", "UHM590.41012" in meta and "UHM590.41012" not in stems)
newst = [s for s in stems if meta[s]["cohort"] != "WF22"]
check("WF23/WF24 stems in tank metadata", all(s in tank for s in newst))
check("egg mass agrees between metadata files", all(em2[s] == meta[s]["egg_mass_cohort"] for s in newst))
check("WF22 design agrees with metadata", all(d22[s][0] == s.split(".")[0] and d22[s][3] == meta[s]["egg_mass_cohort"]
      for s in stems if meta[s]["cohort"] == "WF22"))
# contrasts: name -> {animal: (group, [stems])}
C = collections.OrderedDict()
for nm, emc, tanks in (("WF24_EM4_tank", "WF24_EM4", ("4.1", "4.2")), ("WF23_EM3_tank", "WF23_EM3", ("5.0", "6.0"))):
    C[nm] = {s.split(".")[0]: (tank[s], [s]) for s in newst if em2[s] == emc}
    if set(g for g, _ in C[nm].values()) != set(tanks): raise SystemExit("ABORT: tanks " + nm)
    if len(C[nm]) != sum(1 for s in newst if em2[s] == emc): raise SystemExit("ABORT: animal with two samples " + nm)
w22 = [s for s in stems if meta[s]["cohort"] == "WF22"]
for nm, months in (("WF22_EM_month2", (2,)), ("WF22_EM_month3", (3,)), ("WF22_EM_months23_repro", (2, 3))):
    D = {}
    for s in w22:
        if d22[s][1] in months:
            a = d22[s][0]; D.setdefault(a, (d22[s][3], []))[1].append(s)
    C[nm] = D
    check("%s: 15 animals, %d sample(s) each" % (nm, len(months)), len(D) == 15 and all(len(v[1]) == len(months) for v in D.values()))
print("\nanimals per tank after exclusions (all tanks):")
cnt = collections.defaultdict(list)
for s in newst: cnt[(em2[s], tank[s])].append(s)
for (e, t), L in sorted(cnt.items()):
    tr = collections.Counter(trt2[s] for s in L)
    print("  %-9s tank %s animals %2d treatments %s" % (e, t, len(L), dict(sorted(tr.items()))))
    summ.append(("animals_per_tank", e, "tank " + t, len(L), "treatments " + ";".join("%s:%d" % kv for kv in sorted(tr.items()))))
for e, ts in (("WF24_EM4", ("4.1", "4.2")), ("WF23_EM3", ("5.0", "6.0"))):
    a, b = [set(trt2[s] for s in cnt[(e, t)]) for t in ts]
    check("%s every treatment in both tanks" % e, a == b, "only %s: %s; only %s: %s" % (ts[0], sorted(a - b), ts[1], sorted(b - a)))
for nm in C:
    g = collections.Counter(v[0] for v in C[nm].values())
    print("  contrast %-24s %s" % (nm, dict(sorted(g.items()))))
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
si = {s: j for j, s in enumerate(sorted(stems))}
S = {(nm, m): np.zeros((len(L), len(si)), dtype=bool) for nm, L in SETS.items() for m in ("assembly", "reads>=5")}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t"); x = idx.get(r)
        if x is None: continue
        j = si.get(m.split("__", 1)[0])
        if j is not None: S[(x[0], "assembly")][x[1], j] = True
seen = 0
with gzip.open(R + "count_matrix_primary_n88.tsv.gz", "rt") as fh:
    cols = next(fh).rstrip("\n").split("\t")[1:]
    cj = [(j, si[s]) for j, s in enumerate(cols) if s in si]
    for l in fh:
        r, rest = l.split("\t", 1); x = idx.get(r)
        if x is None: continue
        seen += 1; v = rest.rstrip("\n").split("\t")
        for j, k in cj:
            if float(v[j]) >= 5: S[(x[0], "reads>=5")][x[1], k] = True
check("units found in count matrix", seen == len(idx), "%d of %d" % (seen, len(idx)))
check("all retained stems in count matrix", len(cj) == len(si), "%d of %d" % (len(cj), len(si)))
def restricted(M, lab):
    tot = M.sum(1); out = {}
    for e in sorted(set(lab)):
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1)
        out[e] = (n >= 3) & (n >= 0.5 * g.sum()) & (tot - n == 0)
    return out
rows = []; obs_sets = {}; ratio = {}
for nm, D in C.items():
    animals = sorted(D); lab0 = [D[a][0] for a in animals]; G = sorted(set(lab0))
    rng = np.random.default_rng(34); perms = [list(rng.permutation(lab0)) for _ in range(NPERM)]
    ndist = len(set(tuple(p) for p in perms))
    for mode in ("assembly", "reads>=5"):
        A = {}
        for st in ("EU", "K"):
            M = S[(st, mode)]
            A[st] = np.stack([M[:, [si[s] for s in D[a][1]]].any(1) for a in animals], 1)
        oE, oK = restricted(A["EU"], lab0), restricted(A["K"], lab0)
        for g in G: obs_sets[(nm, mode, g)] = oE[g]
        obs = {g: (int(oE[g].sum()), int(oK[g].sum())) for g in G}; obs["all"] = tuple(sum(obs[g][i] for g in G) for i in (0, 1))
        N = {g: [] for g in G + ["all"]}
        for lab in perms:
            cE, cK = restricted(A["EU"], lab), restricted(A["K"], lab)
            t = {g: (int(cE[g].sum()), int(cK[g].sum())) for g in G}; t["all"] = tuple(sum(t[g][i] for g in G) for i in (0, 1))
            for g in t: N[g].append(t[g])
        print("\n%s, %s (animals %d; distinct labelings among permutations %d)" % (nm, mode, len(animals), ndist))
        for g in G + ["all"]:
            na = len(animals) if g == "all" else lab0.count(g)
            ne = np.array([x[0] for x in N[g]]); nk = np.array([x[1] for x in N[g]]); nd = ne - nk
            vals = (("EU", obs[g][0], ne), ("K", obs[g][1], nk), ("D_EU_minus_K", obs[g][0] - obs[g][1], nd))
            out = []
            for st, o, n in vals:
                p = (1 + int((n >= o).sum())) / (NPERM + 1)
                r = "%.2f" % (o / n.mean()) if st != "D_EU_minus_K" and n.mean() > 0 else "NA"
                rows.append((nm, mode, g, na, st, o, "%.1f" % n.mean(), "%.1f" % n.std(), "%.4f" % p, r, ndist))
                if g == "all" and st != "D_EU_minus_K": ratio[(nm, mode, st)] = (o, n.mean(), r, p)
                out.append("%s %d (null %.1f sd %.1f, p %.4f, ratio %s)" % (st, o, n.mean(), n.std(), p, r))
            print("  %-9s n=%2d  %s" % (g, na, " | ".join(out)))
# reproduction check vs Step 34c3b (assembly, all egg masses)
rr = {(r[4]): r for r in rows if r[0] == "WF22_EM_months23_repro" and r[1] == "assembly" and r[2] == "all"}
check("reproduces 34c3b WF22 months 2-3 (EU 948, K 616, D 332, null D -7.1)",
      rr["EU"][5] == 948 and rr["K"][5] == 616 and rr["D_EU_minus_K"][5] == 332 and rr["D_EU_minus_K"][6] == "-7.1",
      "got EU %s K %s D %s null D %s" % (rr["EU"][5], rr["K"][5], rr["D_EU_minus_K"][5], rr["D_EU_minus_K"][6]))
# overlap
print("\noverlap of tank-restricted EU reps with WF22 egg-mass-restricted EU reps:")
for mode in ("assembly", "reads>=5"):
    for ref in ("WF22_EM_months23_repro", "WF22_EM_month2", "WF22_EM_month3"):
        wr = np.zeros(len(EU), bool)
        for (nm, md, g), v in obs_sets.items():
            if nm == ref and md == mode: wr |= v
        for nm in ("WF24_EM4_tank", "WF23_EM3_tank"):
            tr = np.zeros(len(EU), bool)
            for (n2, md, g), v in obs_sets.items():
                if n2 == nm and md == mode: tr |= v
            ov = int((tr & wr).sum())
            print("  %-9s %-14s restricted %4d vs %-24s restricted %4d: overlap %d" % (mode, nm, int(tr.sum()), ref, int(wr.sum()), ov))
            summ.append(("overlap_EU", mode, "%s vs %s" % (nm, ref), ov, "tank-restricted %d; WF22-restricted %d" % (int(tr.sum()), int(wr.sum()))))
print("\nobserved / null-mean ratios, all groups summed (EU | K):")
for mode in ("assembly", "reads>=5"):
    for nm in C:
        e, k = ratio[(nm, mode, "EU")], ratio[(nm, mode, "K")]
        print("  %-9s %-24s EU %4d / %6.1f = %-5s (p %.4f) | K %4d / %6.1f = %-5s (p %.4f)" % (mode, nm, e[0], e[1], e[2], e[3], k[0], k[1], k[2], k[3]))
        summ.append(("ratio_obs_over_null", mode, nm, "EU %s; K %s" % (e[2], k[2]),
                     "EU %d / %.1f p %.4f; K %d / %.1f p %.4f" % (e[0], e[1], e[3], k[0], k[1], k[3])))
for c in checks: summ.append(("check", c[1], c[0], "", c[2]))
with open(OUT1, "w") as o:
    o.write("contrast\tpresence\tgroup\tn_animals\tset\tobserved_restricted\tnull_mean\tnull_sd\tp_one_sided\tobs_over_null_mean\tdistinct_permutations\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
with open(OUT2, "w") as o:
    o.write("item\tpresence_or_status\tunit\tvalue\tnote\n")
    for r in summ: o.write("\t".join(map(str, r)) + "\n")
print("wrote", OUT1, len(rows), "rows;", OUT2, len(summ), "rows")
print("checks failed:", sum(1 for c in checks if c[1] == "FAIL"))
