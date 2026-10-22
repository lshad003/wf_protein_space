# Step 34t: (1) WF22 months 2-3 egg-mass prevalence breadth and (2) tank restriction at c_id95, for Fig. 3.
# Units as scripts/step34c3_eggmass.py / step34m_tank_restriction.py: 41,114 strict_final EU reps and 41,114 random
#   strict_final K reps (random.seed(343)); animals sorted by ID; presence in an animal = any assembled member
#   (results/clusters_95_n88.tsv) in any retained sample, or >= 5 reads at c_id95 (results/step34r/matrix_c_id95_n88.tsv.gz,
#   scripts/step34r_matrices.py).
# Rule (Step 34c3): a rep meets the rule in group g when present in >= 50% of g's animals AND >= 3 animals; restricted =
#   meets it in g and absent from every other animal of the contrast.
# 1. WF22 months 2-3 pooled (15 animals, egg masses 3/3/9; metadata/wf22_design.tsv), assembly. Reproduction check:
#    restricted EU 948, K 616 (Step 34c3b). Then reps meeting the rule (presence part only, no absence condition) in
#    exactly 1, 2 and 3 egg masses; for exactly 1, how many are also restricted.
# 2. Tanks within one egg mass: WF24 EM4 4.1 vs 4.2, WF23 EM3 5.0 vs 6.0 (metadata/wf23_wf24_sequenced_metadata.tsv;
#    UHM586.41010, UHM590.41012 excluded), assembly and >= 5 reads at c_id95, both tanks summed. Null: tank labels permuted
#    among the egg mass's animals (sizes kept), 1,000 permutations, fresh numpy default_rng(34) per contrast, EU and K on
#    the same permutations; one-sided p = (1 + #null >= obs) / 1,001 for EU, K and D = EU - K.
#    Reproduction checks: WF24 assembly EU 822, K 822, WF23 assembly EU 673, K 1,200 (Step 34m); WF24 c_id95 EU 1,064,
#    K 866, D 198 (Step 34r).
# Aborts before writing if any reproduction check fails.
# Outputs (refuse to overwrite): results/step34t/wf22_eggmass_breadth.tsv, results/step34t/tank_restriction_c_id95.tsv
import os, gzip, random, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; O = R + "step34t/"
OUT1, OUT2 = O + "wf22_eggmass_breadth.tsv", O + "tank_restriction_c_id95.tsv"
for p in (OUT1, OUT2):
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
d22 = {}
for l in list(open(W + "/metadata/wf22_design.tsv"))[1:]:
    s, a, t, e, m = l.rstrip("\n").split("\t"); d22[s] = (a, int(m), "WF22_EM" + e)
stems = sorted(s for s in meta if s not in excl)
C = collections.OrderedDict()
D = {}
for s in stems:
    if meta[s]["cohort"] == "WF22" and d22[s][1] in (2, 3): D.setdefault(d22[s][0], (d22[s][2], []))[1].append(s)
C["WF22_EM_months23"] = D
for nm, emc in (("WF24_EM4_tank", "WF24_EM4"), ("WF23_EM3_tank", "WF23_EM3")):
    C[nm] = {s.split(".")[0]: (tank[s], [s]) for s in stems if meta[s]["cohort"] != "WF22" and em2[s] == emc}
for nm, D in C.items():
    print("%s: animals %d, groups %s" % (nm, len(D), dict(collections.Counter(v[0] for v in D.values()))), flush=True)
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
if seen != len(idx) or len(cj) != len(si): raise SystemExit("ABORT: matrix coverage %d/%d reps, %d/%d stems" % (seen, len(idx), len(cj), len(si)))
def animal_matrix(st, mode, D, animals):
    M = S[(st, mode)]; return np.stack([M[:, [si[s] for s in D[a][1]]].any(1) for a in animals], 1)
def meets(M, lab):
    out = {}
    for e in sorted(set(lab)):
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1); out[e] = (n >= 3) & (n >= 0.5 * g.sum())
    return out
def restricted(M, lab):
    tot = M.sum(1); out = {}
    for e in sorted(set(lab)):
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1)
        out[e] = (n >= 3) & (n >= 0.5 * g.sum()) & (tot - n == 0)
    return out
bad = []
# 1
D = C["WF22_EM_months23"]; an = sorted(D); lab = [D[a][0] for a in an]
rows1 = []
for st in ("EU", "K"):
    A = animal_matrix(st, "assembly", D, an); rs = restricted(A, lab); mt = meets(A, lab)
    nr = sum(int(v.sum()) for v in rs.values()); nm_ = np.sum([mt[e] for e in mt], 0)
    anyr = np.any([rs[e] for e in rs], 0)
    exp = {"EU": 948, "K": 616}[st]
    print("1. %s restricted %d (expected %d)" % (st, nr, exp), flush=True)
    if nr != exp: bad.append("WF22 restricted %s %d != %d" % (st, nr, exp))
    for k in (1, 2, 3):
        sel = nm_ == k
        rows1.append((st, k, int(sel.sum()), int((sel & anyr).sum()), len(A), nr))
        print("   %s meets rule in exactly %d egg mass(es): %d (also restricted %d)" % (st, k, int(sel.sum()), int((sel & anyr).sum())))
    for e in sorted(mt): rows1.append((st, "only_" + e, int((mt[e] & (nm_ == 1)).sum()), int(rs[e].sum()), len(A), nr))
# 2
rows2 = []
for nm in ("WF24_EM4_tank", "WF23_EM3_tank"):
    D = C[nm]; an = sorted(D); lab0 = [D[a][0] for a in an]; G = sorted(set(lab0))
    rng = np.random.default_rng(34); perms = [list(rng.permutation(lab0)) for _ in range(NPERM)]
    ndist = len(set(tuple(p) for p in perms))
    for mode in ("assembly", "reads>=5_c_id95"):
        AE, AK = animal_matrix("EU", mode, D, an), animal_matrix("K", mode, D, an)
        sc = lambda lab: (sum(int(v.sum()) for v in restricted(AE, lab).values()), sum(int(v.sum()) for v in restricted(AK, lab).values()))
        oE, oK = sc(lab0); N = np.array([sc(p) for p in perms]); nd = N[:, 0] - N[:, 1]
        for st, o, n in (("EU", oE, N[:, 0]), ("K", oK, N[:, 1]), ("D_EU_minus_K", oE - oK, nd)):
            p = (1 + int((n >= o).sum())) / (NPERM + 1)
            rows2.append((nm, mode, st, len(an), ";".join("%s:%d" % (g, lab0.count(g)) for g in G), o, "%.1f" % n.mean(), "%.1f" % n.std(), "%.4f" % p, ndist))
            print("2. %-14s %-16s %-13s obs %5d null %7.1f +/- %6.1f p %.4f (distinct labelings %d)" % (nm, mode, st, o, n.mean(), n.std(), p, ndist), flush=True)
chk = {(r[0], r[1], r[2]): r[5] for r in rows2}
for key, exp in ((("WF24_EM4_tank", "assembly", "EU"), 822), (("WF24_EM4_tank", "assembly", "K"), 822),
                 (("WF23_EM3_tank", "assembly", "EU"), 673), (("WF23_EM3_tank", "assembly", "K"), 1200),
                 (("WF24_EM4_tank", "reads>=5_c_id95", "EU"), 1064), (("WF24_EM4_tank", "reads>=5_c_id95", "K"), 866),
                 (("WF24_EM4_tank", "reads>=5_c_id95", "D_EU_minus_K"), 198)):
    print("CHECK %-55s got %d expected %d %s" % (" ".join(key), chk[key], exp, "PASS" if chk[key] == exp else "FAIL"))
    if chk[key] != exp: bad.append("%s %d != %d" % (" ".join(key), chk[key], exp))
if bad: raise SystemExit("ABORT, reproduction failed, nothing written: %s" % bad)
os.makedirs(O, exist_ok=True)
with open(OUT1, "w") as o:
    o.write("set\tn_egg_masses_meeting_rule\tn_reps\tof_which_restricted\tn_reps_tested\trestricted_total\n")
    for r in rows1: o.write("\t".join(map(str, r)) + "\n")
with open(OUT2, "w") as o:
    o.write("contrast\tpresence\tset\tn_animals\tanimals_per_tank\tobserved_restricted\tnull_mean\tnull_sd\tp_one_sided\tdistinct_permutations\n")
    for r in rows2: o.write("\t".join(map(str, r)) + "\n")
print("wrote", OUT1, OUT2)
