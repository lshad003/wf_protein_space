# Step 34c3b: two checks on Step 34c3 (no new data; same units, presence rules and WF22 months 2-3 as step34c3_eggmass.py).
# 1. Paired permutation test, WF22: statistic D = (EU restricted) - (K restricted), summed over egg masses and per egg
#    mass. 1,000 within-year permutations of egg-mass labels among WF22 animals (group sizes kept, seed 34), each
#    permutation applied to EU and K together. Step 34c3 drew separate permutations per set; this run uses one shared
#    set, so the separate EU and K null means are reported again from it. Empirical p = (1 + #null >= obs) / 1,001.
# 2. EM3 treatment confound: for the EM3-restricted EU reps (837) and K reps (576), the fraction present in each EM3
#    animal, and the mean by treatment (STP1710.7, control, STP1717.1).
# Both presence modes (assembly primary, reads >= 5 sensitivity). Output: results/step34c3b_checks.tsv; refuses to overwrite.
import os, gzip, random, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = R + "step34c3b_checks.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
NPERM = 1000
meta = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); meta[f["stem"]] = f
d22 = {}
with open(W + "/metadata/wf22_design.tsv") as fh:
    next(fh)
    for l in fh:
        s, a, t, e, m = l.rstrip("\n").split("\t"); d22[s] = (a, int(m), t)
keep = {s for s in meta if meta[s]["cohort"] == "WF22" and d22[s][1] in (2, 3)}
animal = {s: s.split(".", 1)[0] for s in keep}
em = {}; trt = {}
for s in keep:
    a = animal[s]; em[a] = meta[s]["egg_mass_cohort"]
    if trt.get(a, d22[s][2]) != d22[s][2]: raise SystemExit("ABORT: animal with two treatments " + a)
    trt[a] = d22[s][2]
animals = sorted(em); ai = {a: i for i, a in enumerate(animals)}
print("WF22 animals %d, samples %d (months 2-3)" % (len(animals), len(keep)))
for a in animals: print("  %-7s %s %s" % (a, em[a], trt[a]))
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
P = {(nm, m): np.zeros((len(L), len(animals)), dtype=bool) for nm, L in SETS.items() for m in ("assembly", "reads>=5")}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t"); x = idx.get(r)
        if x is None: continue
        s = m.split("__", 1)[0]
        if s in keep: P[(x[0], "assembly")][x[1], ai[animal[s]]] = True
with gzip.open(R + "count_matrix_primary_n88.tsv.gz", "rt") as fh:
    cols = next(fh).rstrip("\n").split("\t")[1:]
    cj = [(j, ai[animal[s]]) for j, s in enumerate(cols) if s in keep]
    for l in fh:
        r, rest = l.split("\t", 1); x = idx.get(r)
        if x is None: continue
        v = rest.rstrip("\n").split("\t")
        for j, a in cj:
            if float(v[j]) >= 5: P[(x[0], "reads>=5")][x[1], a] = True
def restricted(M, lab):
    ems = sorted(set(lab)); tot = M.sum(1); out = {}
    for e in ems:
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1); na = int(g.sum())
        out[e] = (n >= 3) & (n >= 0.5 * na) & (tot - n == 0)
    return out
lab0 = [em[a] for a in animals]; EMS = sorted(set(lab0))
rng = np.random.default_rng(34); perms = [list(rng.permutation(lab0)) for _ in range(NPERM)]
rows = []
for mode in ("assembly", "reads>=5"):
    oE, oK = restricted(P[("EU", mode)], lab0), restricted(P[("K", mode)], lab0)
    obs = {e: (int(oE[e].sum()), int(oK[e].sum())) for e in EMS}; obs["all"] = tuple(sum(obs[e][i] for e in EMS) for i in (0, 1))
    N = {e: [] for e in EMS + ["all"]}
    for lab in perms:
        cE, cK = restricted(P[("EU", mode)], lab), restricted(P[("K", mode)], lab)
        t = {e: (int(cE[e].sum()), int(cK[e].sum())) for e in EMS}; t["all"] = tuple(sum(t[e][i] for e in EMS) for i in (0, 1))
        for e in t: N[e].append(t[e])
    print("\n1. paired test WF22, %s: D = EU restricted - K restricted" % mode)
    for e in EMS + ["all"]:
        d0 = obs[e][0] - obs[e][1]; nd = np.array([x[0] - x[1] for x in N[e]])
        ne = np.array([x[0] for x in N[e]]); nk = np.array([x[1] for x in N[e]])
        p = (1 + int((nd >= d0).sum())) / (NPERM + 1); p2 = (1 + int((np.abs(nd - nd.mean()) >= abs(d0 - nd.mean())).sum())) / (NPERM + 1)
        print("   %-9s EU %5d (null %.1f) K %5d (null %.1f) | D obs %5d, null mean %.1f sd %.1f, p (D >= obs) %.4f, two-sided %.4f"
              % (e, obs[e][0], ne.mean(), obs[e][1], nk.mean(), d0, nd.mean(), nd.std(), p, p2))
        rows.append(("paired_test", mode, e, "EU_restricted", obs[e][0], "%.1f" % ne.mean(), ""))
        rows.append(("paired_test", mode, e, "K_restricted", obs[e][1], "%.1f" % nk.mean(), ""))
        rows.append(("paired_test", mode, e, "D_EU_minus_K", d0, "%.1f (sd %.1f)" % (nd.mean(), nd.std()), "p_one_sided %.4f; p_two_sided %.4f" % (p, p2)))
    # 2. EM3 treatment confound
    A3 = [a for a in animals if em[a] == "WF22_EM3"]
    print("\n2. EM3-restricted reps present per EM3 animal, %s" % mode)
    for nm, o in (("EU", oE), ("K", oK)):
        sel = o["WF22_EM3"]; M = P[(nm, mode)][sel]; n = int(sel.sum())
        print("   %s: %d EM3-restricted reps" % (nm, n))
        fr = {}
        for a in A3:
            fr[a] = M[:, ai[a]].mean() if n else float("nan")
            print("     %-7s %-10s %.3f" % (a, trt[a], fr[a]))
            rows.append(("em3_animal", mode, a, nm + " " + trt[a], n, "%.4f" % fr[a], ""))
        for t in ("STP1710.7", "control", "STP1717.1"):
            v = [fr[a] for a in A3 if trt[a] == t]
            print("     %-10s animals %d mean fraction %.3f" % (t, len(v), np.mean(v)))
            rows.append(("em3_treatment", mode, t, nm, len(v), "%.4f" % np.mean(v), "min %.3f max %.3f" % (min(v), max(v))))
        cnt = M[:, [ai[a] for a in A3]].sum(1)
        only = {t: int((M[:, [ai[a] for a in A3 if trt[a] != t]].sum(1) == 0).sum()) for t in ("STP1710.7", "control", "STP1717.1")}
        print("     reps present in all %d EM3 animals: %d | absent from every animal outside one treatment group: %s" % (len(A3), int((cnt == len(A3)).sum()), only))
        rows.append(("em3_summary", mode, "all_EM3_animals", nm, n, int((cnt == len(A3)).sum()), "confined to one treatment: " + ";".join("%s %d" % kv for kv in only.items())))
with open(OUT, "w") as o:
    o.write("check\tpresence\tunit\tmeasure\tn\tvalue\tnote\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("wrote", OUT, "rows:", len(rows))
