# Step 34t5: per-class, per-egg-mass nulls for Fig. 3d (WF22 egg-mass restriction, assembly presence).
# Exactly Step 34c3 (scripts/step34c3_eggmass.py): WF22 months shared by all animals (2, 3), 15 froglets (EM1 3, EM2 3,
#   EM3 9); 41,114 strict_final EU and 41,114 random K reps (random.seed(343)); presence = any assembled member
#   (results/clusters_95_n88.tsv) in a retained sample; restricted = >= 50% of the egg mass's froglets and >= 3, absent
#   from every froglet of the other WF22 egg masses.
# Null "step34c3": egg-mass labels permuted among the WF22 froglets (group sizes kept), 1,000 permutations per set from
#   one numpy default_rng(34) stream advanced in Step 34c3's order (EU assembly, EU reads >= 5, K assembly, ...), so the
#   EU assembly null is permutations 1-1,000 and the K assembly null 2,001-3,000. p = (1 + #null >= obs) / 1,001.
# Null "shared" (Step 34c3b): one default_rng(34) set of 1,000 permutations applied to EU and K together; gives the
#   EU - K test (D 332). Reported alongside for comparison.
# 95th percentile: numpy linear interpolation.
# Reproduction (abort before writing): restricted EU 948 = 48/63/837, K 616 = 29/11/576; Step 34c3 null means
#   EU 4.9/5.8/56.8, K 4.2/3.9/69.6 and p EU 0.0120/0.0140/0.0010, K 0.0310/0.0799/0.0010; D 332, null -7.1, p 0.0010.
# Output (refuses to overwrite): results/step34t/wf22_eggmass_nulls.tsv
import os, random, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; OUT = R + "step34t/wf22_eggmass_nulls.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
NPERM = 1000
meta = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); meta[f["stem"]] = f
d22 = {}
for l in list(open(W + "/metadata/wf22_design.tsv"))[1:]:
    s, a, t, e, m = l.rstrip("\n").split("\t"); d22[s] = (a, int(m), "WF22_EM" + e)
w22 = [s for s in meta if meta[s]["cohort"] == "WF22"]
mon = collections.defaultdict(set)
for s in w22: mon[d22[s][0]].add(d22[s][1])
shared = set.intersection(*mon.values())
if shared != {2, 3}: raise SystemExit("ABORT: shared months %s" % shared)
keep = [s for s in w22 if d22[s][1] in shared]
an = sorted({d22[s][0] for s in keep}); ai = {a: i for i, a in enumerate(an)}
lab0 = [next(d22[s][2] for s in keep if d22[s][0] == a) for a in an]
EMS = sorted(set(lab0))
print("WF22 froglets %d, samples %d, egg masses %s" % (len(an), len(keep), dict(collections.Counter(lab0))))
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
ks = set(keep); st2a = {s: ai[d22[s][0]] for s in keep}
P = {nm: np.zeros((len(L), len(an)), dtype=bool) for nm, L in SETS.items()}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t"); x = idx.get(r)
        if x is None: continue
        s = m.split("__", 1)[0]
        if s in ks: P[x[0]][x[1], st2a[s]] = True
def restricted(M, lab):
    tot = M.sum(1); out = {}
    for e in EMS:
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1); out[e] = int(((n >= 3) & (n >= 0.5 * g.sum()) & (tot - n == 0)).sum())
    out["all"] = sum(out[e] for e in EMS); return out
obs = {st: restricted(P[st], lab0) for st in SETS}
rng = np.random.default_rng(34); seq = [list(rng.permutation(lab0)) for _ in range(4 * NPERM)]
PERMS = {("step34c3", "EU"): seq[0:NPERM], ("step34c3", "K"): seq[2 * NPERM:3 * NPERM]}
rng = np.random.default_rng(34); sh = [list(rng.permutation(lab0)) for _ in range(NPERM)]
PERMS[("shared", "EU")] = sh; PERMS[("shared", "K")] = sh
NULL = {k: [restricted(P[k[1]], p) for p in v] for k, v in PERMS.items()}
rows = []
for meth in ("step34c3", "shared"):
    for st in ("EU", "K"):
        for e in EMS + ["all"]:
            n = np.array([x[e] for x in NULL[(meth, st)]]); o = obs[st][e]
            rows.append((meth, st, e, lab0.count(e) if e != "all" else len(lab0), o, "%.2f" % n.mean(), "%.2f" % n.std(), "%.2f" % np.percentile(n, 95), "%.4f" % ((1 + int((n >= o).sum())) / (NPERM + 1))))
    if meth == "shared":
        for e in EMS + ["all"]:
            n = np.array([x[e] for x in NULL[(meth, "EU")]]) - np.array([x[e] for x in NULL[(meth, "K")]]); o = obs["EU"][e] - obs["K"][e]
            rows.append((meth, "D_EU_minus_K", e, lab0.count(e) if e != "all" else len(lab0), o, "%.2f" % n.mean(), "%.2f" % n.std(), "%.2f" % np.percentile(n, 95), "%.4f" % ((1 + int((n >= o).sum())) / (NPERM + 1))))
print("\n%-9s %-13s %-9s %3s %5s %9s %8s %8s %7s" % ("null", "set", "egg_mass", "n", "obs", "mean", "sd", "p95", "p"))
for r in rows: print("%-9s %-13s %-9s %3s %5s %9s %8s %8s %7s" % r)
T = {(r[0], r[1], r[2]): r for r in rows}; bad = []
def ck(name, got, exp):
    print("CHECK %-48s got %-8s expected %-8s %s" % (name, got, exp, "PASS" if str(got) == str(exp) else "FAIL"))
    if str(got) != str(exp): bad.append(name)
for st, ex in (("EU", (48, 63, 837, 948)), ("K", (29, 11, 576, 616))):
    for e, v in zip(EMS + ["all"], ex): ck("restricted %s %s" % (st, e), obs[st][e], v)
for st, means, ps in (("EU", ("4.9", "5.8", "56.8"), ("0.0120", "0.0140", "0.0010")), ("K", ("4.2", "3.9", "69.6"), ("0.0310", "0.0799", "0.0010"))):
    for e, m, p in zip(EMS, means, ps):
        ck("34c3 null mean %s %s" % (st, e), "%.1f" % np.mean([x[e] for x in NULL[("step34c3", st)]]), m); ck("34c3 p %s %s" % (st, e), T[("step34c3", st, e)][8], p)
ck("D all", T[("shared", "D_EU_minus_K", "all")][4], 332); dn = np.array([x["all"] for x in NULL[("shared", "EU")]]) - np.array([x["all"] for x in NULL[("shared", "K")]])
ck("D null mean (raw mean, one decimal, as 34c3b)", "%.1f" % dn.mean(), "-7.1")
ck("D p", T[("shared", "D_EU_minus_K", "all")][8], "0.0010")
if bad: raise SystemExit("ABORT, reproduction failed, nothing written: %s" % bad)
with open(OUT, "w") as o:
    o.write("null\tset\tegg_mass\tn_froglets\tobserved_restricted\tnull_mean\tnull_sd\tnull_p95\tp_one_sided\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("wrote", OUT)
