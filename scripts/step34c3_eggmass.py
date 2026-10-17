# Step 34c3: egg-mass-level recurrence of strict_final EU reps (descriptive, with a permutation null).
# Metadata: metadata/samples_88_treatment.tsv (stem, cohort, treatment, egg_mass, egg_mass_cohort); WF22 animal and
#   month from metadata/wf22_design.tsv; WF23/WF24 one sample per animal (animal = stem prefix before ".").
#   WF22 uses only the months every WF22 animal shares (intersection over animals).
# Units: strict_final EU reps (results/step34a_rep_origin_v2_final.tsv, strict_final = 1; 41,114); baseline an equal
#   number of strict_final K reps drawn at random (seed 343).
# Presence in an animal: primary = any member protein (results/clusters_95_n88.tsv) assembled from that animal's
#   retained samples (as Step 34c); sensitivity = >= 5 reads in any retained sample (results/count_matrix_primary_n88.tsv.gz).
# Within-egg-mass recurrent (rep, egg mass): present in >= 50% of the egg mass's animals AND >= 3 animals.
# Egg-mass-restricted: within-egg-mass recurrent AND absent from every animal of the other egg masses that year.
# Null: egg-mass labels permuted among animals within the year, group sizes kept, 1,000 permutations (seed 34); per egg
#   mass and per year, counts of within-recurrent and restricted reps; empirical p = (1 + #null >= observed) / 1,001.
# WF22 primary; WF24 descriptive (EM2 has 3 animals); WF23 one egg mass, no contrast (within-recurrent only).
# Confounds: STP1710.7 occurs only in WF22 EM3; in WF24 egg mass equals tank.
# Overlap: within-egg-mass recurrent EU reps vs reps of the 4,509 strict_final cross-year recurrent families
#   (results/step34c2c_strict_v2_families.tsv recurrent = 1, excluded = 0; family reps from step34c_strict_eu_clusters30.tsv).
# Output: results/step34c3_eggmass_recurrence.tsv; refuses to overwrite.
import os, gzip, random, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = R + "step34c3_eggmass_recurrence.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
NPERM, SEED_K, SEED_P = 1000, 343, 34
meta = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); meta[f["stem"]] = f
if len(meta) != 88: raise SystemExit("ABORT: %d stems" % len(meta))
d22 = {}
with open(W + "/metadata/wf22_design.tsv") as fh:
    next(fh)
    for l in fh:
        s, a, t, e, m = l.rstrip("\n").split("\t"); d22[s] = (a, int(m), e)
w22 = [s for s in meta if meta[s]["cohort"] == "WF22"]
if set(w22) - set(d22): raise SystemExit("ABORT: WF22 stems missing from wf22_design.tsv")
bad = [s for s in w22 if d22[s][0] != s.split(".")[0] or "WF22_EM" + d22[s][2] != meta[s]["egg_mass_cohort"]]
if bad: raise SystemExit("ABORT: WF22 design disagrees with metadata: %s" % bad)
mon = collections.defaultdict(set)
for s in w22: mon[d22[s][0]].add(d22[s][1])
shared = set.intersection(*mon.values())
print("WF22 months per animal:", {a: sorted(m) for a, m in sorted(mon.items())})
print("WF22 months shared by all %d animals: %s" % (len(mon), sorted(shared)))
keep = {s for s in meta if meta[s]["cohort"] != "WF22" or d22[s][1] in shared}
animal = {s: s.split(".", 1)[0] for s in keep}
em = {}; year = {}
for s in keep:
    a = animal[s]; e = meta[s]["egg_mass_cohort"]; y = meta[s]["cohort"]
    if em.get(a, e) != e: raise SystemExit("ABORT: animal in two egg masses " + a)
    em[a] = e; year[a] = y
print("\nanimals and retained samples per egg_mass_cohort (WF22 restricted to months %s):" % sorted(shared))
for e in sorted(set(em.values())):
    A = [a for a in em if em[a] == e]; Sx = [s for s in keep if em[animal[s]] == e]
    tr = collections.Counter(meta[s]["treatment"] for s in Sx)
    print("  %-9s animals %2d | samples %2d | treatments (samples) %s" % (e, len(A), len(Sx), dict(tr)))
print("  samples retained %d of 88" % len(keep), flush=True)
# rep sets
EU, K = [], []
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, ifn = h.index("class4"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ifn] == "1":
            if f[ic] == "EU": EU.append(f[0])
            elif f[ic] == "K": K.append(f[0])
print("\nstrict_final EU %d | K %d" % (len(EU), len(K)), flush=True)
if len(EU) != 41114: raise SystemExit("ABORT: expected 41,114 EU")
random.seed(SEED_K); Kb = sorted(random.sample(sorted(K), len(EU)))
SETS = {"EU": sorted(EU), "K": Kb}
idx = {}
for nm, L in SETS.items():
    for i, r in enumerate(L): idx[r] = (nm, i)
animals = sorted(em); ai = {a: i for i, a in enumerate(animals)}
P = {(nm, m): np.zeros((len(L), len(animals)), dtype=bool) for nm, L in SETS.items() for m in ("assembly", "reads>=5")}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        x = idx.get(r)
        if x is None: continue
        s = m.split("__", 1)[0]
        if s in keep: P[(x[0], "assembly")][x[1], ai[animal[s]]] = True
with gzip.open(R + "count_matrix_primary_n88.tsv.gz", "rt") as fh:
    cols = next(fh).rstrip("\n").split("\t")[1:]
    cj = [(j, ai[animal[s]]) for j, s in enumerate(cols) if s in keep]
    seen = 0
    for l in fh:
        r, rest = l.split("\t", 1)
        x = idx.get(r)
        if x is None: continue
        seen += 1; v = rest.rstrip("\n").split("\t")
        for j, a in cj:
            if float(v[j]) >= 5: P[(x[0], "reads>=5")][x[1], a] = True
print("count-matrix rows found for the units: %d of %d" % (seen, len(idx)), flush=True)
if seen != len(idx): raise SystemExit("ABORT: units missing from count matrix")
for k, M in P.items():
    print("  %s %-9s reps present in >= 1 retained animal: %d" % (k[0], k[1], int(M.any(1).sum())))
# recurrence
def counts(M, lab, yr):
    A = [i for i, a in enumerate(animals) if year[a] == yr]
    sub = M[:, A]; L = [lab[i] for i in A]; ems = sorted(set(L))
    G = {e: np.array([x == e for x in L]) for e in ems}
    n = {e: sub[:, G[e]].sum(1) for e in ems}
    out = {}; tot = sub.sum(1)
    for e in ems:
        na = int(G[e].sum()); w = (n[e] >= 3) & (n[e] >= 0.5 * na)
        rs = w & (tot - n[e] == 0)
        out[e] = (w, rs)
    return out
rows = []; rng = np.random.default_rng(SEED_P)
within_eu = {}
for yr in ("WF22", "WF24", "WF23"):
    A = [a for a in animals if year[a] == yr]; lab0 = [em[a] for a in animals]
    for (nm, mode), M in sorted(P.items()):
        obs = counts(M, lab0, yr)
        if nm == "EU" and mode == "assembly":
            for e, (w, rs) in obs.items(): within_eu[e] = (w, rs)
        ems = sorted(obs)
        nw = {e: int(obs[e][0].sum()) for e in ems}; nr = {e: int(obs[e][1].sum()) for e in ems}
        if len(ems) > 1:
            Ai = [ai[a] for a in A]; labs = [em[a] for a in A]
            NW = {e: [] for e in ems + ["all"]}; NR = {e: [] for e in ems + ["all"]}
            for _ in range(NPERM):
                perm = list(rng.permutation(labs)); lab = list(lab0)
                for i, x in zip(Ai, perm): lab[i] = x
                c = counts(M, lab, yr)
                for e in ems: NW[e].append(int(c[e][0].sum())); NR[e].append(int(c[e][1].sum()))
                NW["all"].append(sum(int(c[e][0].sum()) for e in ems)); NR["all"].append(sum(int(c[e][1].sum()) for e in ems))
            nw["all"] = sum(nw[e] for e in ems); nr["all"] = sum(nr[e] for e in ems)
            for e in ems + ["all"]:
                pw = (1 + sum(x >= nw[e] for x in NW[e])) / (NPERM + 1); pr = (1 + sum(x >= nr[e] for x in NR[e])) / (NPERM + 1)
                na = sum(1 for a in A if e == "all" or em[a] == e)
                rows.append((yr, e, nm, mode, na, nw[e], "%.1f" % np.mean(NW[e]), "%.4f" % pw, nr[e], "%.1f" % np.mean(NR[e]), "%.1f" % np.std(NR[e]), "%.4f" % pr))
        else:
            for e in ems:
                rows.append((yr, e, nm, mode, len(A), nw[e], "NA", "NA", "NA", "NA", "NA", "NA"))
        print("done %s %s %s" % (yr, nm, mode), flush=True)
HDR = ("year", "egg_mass", "set", "presence", "n_animals", "within_recurrent", "null_mean_within", "p_within",
       "restricted", "null_mean_restricted", "null_sd_restricted", "p_restricted")
with open(OUT, "w") as o:
    o.write("\t".join(HDR) + "\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("\n" + "\t".join(HDR))
for r in rows: print("\t".join(map(str, r)))
# overlap with cross-year recurrent families
fam_of = {}
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t"); fam_of[b] = a
rec = set()
with open(R + "step34c2c_strict_v2_families.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ir, ix = h.index("recurrent"), h.index("excluded")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ir] == "1" and f[ix] == "0": rec.add(f[0])
print("\ncross-year recurrent families (strict_final): %d" % len(rec))
if len(rec) != 4509: raise SystemExit("ABORT: expected 4,509")
EUl = SETS["EU"]; inrec = np.array([fam_of[r] in rec for r in EUl])
print("strict_final EU reps in cross-year recurrent families: %d" % int(inrec.sum()))
anyw = np.zeros(len(EUl), bool); anyr = np.zeros(len(EUl), bool)
for e, (w, rs) in sorted(within_eu.items()):
    anyw |= w
    if e.startswith("WF23"): continue
    anyr |= rs
    print("  %-9s within-recurrent %5d, of which in cross-year recurrent families %5d | restricted %5d, of which %d"
          % (e, int(w.sum()), int((w & inrec).sum()), int(rs.sum()), int((rs & inrec).sum())))
print("  any egg mass: within-recurrent %d, in cross-year recurrent families %d (%.1f%%) | restricted (WF22/WF24) %d, in cross-year recurrent %d"
      % (int(anyw.sum()), int((anyw & inrec).sum()), 100.0 * (anyw & inrec).sum() / max(1, anyw.sum()), int(anyr.sum()), int((anyr & inrec).sum())))
print("  cross-year recurrent reps that are within-recurrent in >= 1 egg mass: %d of %d" % (int((anyw & inrec).sum()), int(inrec.sum())))
print("wrote", OUT, "rows:", len(rows))
