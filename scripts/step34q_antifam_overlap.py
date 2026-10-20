# Step 34q: overlap of the Step 34n AntiFam hits (results/step34n_antifam_eu.tblout, --cut_ga; 360 EU reps) with the
#   main EU sets. Restricted lists are rebuilt and saved (they were not saved by Steps 34c3 and 34m):
#   - WF22 egg-mass-restricted (Step 34c3 rules, scripts/step34c3_eggmass.py: WF22 months 2-3, assembly presence, >= 3
#     animals and >= 50% of the egg mass, absent from the other WF22 egg masses); expected EU 948, K 616.
#   - WF24 EM4 tank-restricted (Step 34m rules, scripts/step34m_tank_restriction.py: UHM586/UHM590 excluded, assembly
#     presence, tank 4.1 vs 4.2); expected EU 822 (304 + 518), K 822 (242 + 580).
#   K = the 41,114 random strict_final K reps (seed 343) of Steps 34c3/34m.
# Other sets: the 4,509 recurrent family reps and all 5,107 EU reps of recurrent families (results/step34c2c_strict_v2_
#   families.tsv recurrent 1, excluded 0; results/step34c_strict_eu_clusters30.tsv); the 186 structure-hit families (best
#   E <= 1e-3 in results/step34d_prostt5/queries_vs_{PDB,CATH50}.tsv) and all EU reps of those families; the 30,070
#   NMPFams2 queries (strict_final EU reps in results/unk_filtered_aplakidou.faa; Step 34b2).
# Outputs: results/step34q_antifam_overlap.tsv (per set), results/step34q_antifam_hits_flags.tsv (the 360 hits with a
#   flag per set), results/step34q/{eggmass,tank}_restricted_{EU,K}.tsv (rep, group). Refuses to overwrite.
import os, random, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT1, OUT2 = R + "step34q_antifam_overlap.tsv", R + "step34q_antifam_hits_flags.tsv"
LD = R + "step34q/"
LISTS = [LD + "%s_restricted_%s.tsv" % (a, b) for a in ("eggmass", "tank") for b in ("EU", "K")]
for p in [OUT1, OUT2] + LISTS:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
os.makedirs(LD, exist_ok=True)
meta = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); meta[f["stem"]] = f
month = {}
with open(W + "/metadata/wf22_design.tsv") as fh:
    next(fh)
    for l in fh:
        s, a, t, e, m = l.rstrip("\n").split("\t"); month[s] = int(m)
excl = {l.split("\t")[0] for l in list(open(W + "/metadata/wf24_excluded.tsv"))[1:]}
tank = {}; em2 = {}
with open(W + "/metadata/wf23_wf24_sequenced_metadata.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); tank[f["stem"]] = f["tank"]; em2[f["stem"]] = f["wf_year"] + "_EM" + f["egg_mass"]
# contrasts: stem -> (animal, group)
C = {"eggmass": {s: (s.split(".")[0], meta[s]["egg_mass_cohort"]) for s in meta if meta[s]["cohort"] == "WF22" and month[s] in (2, 3)},
     "tank": {s: (s.split(".")[0], tank[s]) for s in meta if s not in excl and em2.get(s) == "WF24_EM4"}}
EU, K = [], []
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, isf = h.index("class4"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[isf] == "1":
            if f[ic] == "EU": EU.append(f[0])
            elif f[ic] == "K": K.append(f[0])
if len(EU) != 41114: raise SystemExit("ABORT: %d EU" % len(EU))
random.seed(343); Kb = sorted(random.sample(sorted(K), len(EU)))
SETS = {"EU": sorted(EU), "K": Kb}
idx = {r: (nm, i) for nm, L in SETS.items() for i, r in enumerate(L)}
A = {}
for cn, D in C.items():
    animals = sorted({a for a, g in D.values()}); A[cn] = (animals, [next(g for s, (a2, g) in D.items() if a2 == a) for a in animals])
P = {(cn, nm): np.zeros((len(SETS[nm]), len(A[cn][0])), bool) for cn in C for nm in SETS}
ai = {cn: {a: i for i, a in enumerate(A[cn][0])} for cn in C}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t"); x = idx.get(r)
        if x is None: continue
        s = m.split("__", 1)[0]
        for cn, D in C.items():
            if s in D: P[(cn, x[0])][x[1], ai[cn][D[s][0]]] = True
restr = {}
for cn in C:
    lab = A[cn][1]
    for nm in SETS:
        M = P[(cn, nm)]; tot = M.sum(1); out = []
        for g in sorted(set(lab)):
            gm = np.array([x == g for x in lab]); n = M[:, gm].sum(1)
            rs = (n >= 3) & (n >= 0.5 * gm.sum()) & (tot - n == 0)
            out += [(SETS[nm][i], g) for i in np.where(rs)[0]]
        restr[(cn, nm)] = out
        with open(LD + "%s_restricted_%s.tsv" % (cn, nm), "w") as o:
            o.write("rep\tgroup\n")
            for r, g in out: o.write("%s\t%s\n" % (r, g))
        print("%s %s restricted %d by group %s" % (cn, nm, len(out), dict(collections.Counter(g for _, g in out))))
exp = {("eggmass", "EU"): 948, ("eggmass", "K"): 616, ("tank", "EU"): 822, ("tank", "K"): 822}
bad = {k: len(restr[k]) for k in exp if len(restr[k]) != exp[k]}
if bad: raise SystemExit("ABORT: restricted counts differ from Steps 34c3/34m: %s" % bad)
# other sets
fam_of = {}
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t"); fam_of[b] = a
rec = set()
for l in list(open(R + "step34c2c_strict_v2_families.tsv"))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[8] == "1" and f[9] == "0": rec.add(f[0])
eus = set(EU)
rec_reps = {r for r in eus if fam_of.get(r) in rec}
q = [l[1:].split()[0] for l in open(R + "step34d_queries.faa") if l[0] == ">"]
best = {}
for db in ("PDB", "CATH50"):
    for l in open(R + "step34d_prostt5/queries_vs_%s.tsv" % db):
        p = l.split("\t"); e = float(p[2]); best[p[0]] = min(e, best.get(p[0], 9e9))
shit = {x for x in q if best.get(x, 9e9) <= 1e-3}
shit_fams = {fam_of.get(x, x) for x in shit}
shit_reps = {r for r in eus if fam_of.get(r) in shit_fams}
nmp = {l[1:].split()[0] for l in open(R + "unk_filtered_aplakidou.faa") if l[0] == ">"} & eus
hits = set()
for l in open(R + "step34n_antifam_eu.tblout"):
    if l[0] != "#": hits.add(l.split()[0])
hitfam = collections.defaultdict(set)
for l in open(R + "step34n_antifam_eu.tblout"):
    if l[0] != "#": f = l.split(); hitfam[f[0]].add(f[2])
checks = {"recurrent family reps": (len(rec), 4509), "EU reps in recurrent families": (len(rec_reps), 5107),
          "structure-hit families": (len(shit), 186), "NMPFams2 queries": (len(nmp), 30070), "AntiFam hits": (len(hits), 360)}
for k, (a, b) in checks.items():
    print("check %-32s %d (expected %d) %s" % (k, a, b, "PASS" if a == b else "FAIL"))
if any(a != b for a, b in checks.values()): raise SystemExit("ABORT: a set size differs")
S = [("all strict_final EU reps", eus, "Step 34n query set"),
     ("recurrent family reps (4,509)", rec & eus, "family representative of each recurrent family"),
     ("EU reps in recurrent families", rec_reps, "all strict_final EU reps whose 30% family is recurrent"),
     ("WF22 egg-mass-restricted EU (948)", {r for r, g in restr[("eggmass", "EU")]}, "Step 34c3 rules, rebuilt"),
     ("structure-hit families (186)", shit, "family rep with best E <= 1e-3, Step 34d"),
     ("EU reps in structure-hit families", shit_reps, "all EU reps of the 186 families"),
     ("WF24 EM4 tank-restricted EU (822)", {r for r, g in restr[("tank", "EU")]}, "Step 34m rules, rebuilt"),
     ("NMPFams2 queries (30,070)", nmp, "strict_final EU reps in unk_filtered_aplakidou.faa")]
with open(OUT1, "w") as o:
    o.write("set\tset_size\tantifam_hits_in_set\tpct_of_set\tpct_of_360_hits\tnote\n")
    for nm, st, note in S:
        k = len(hits & st)
        o.write("%s\t%d\t%d\t%.2f\t%.2f\t%s\n" % (nm, len(st), k, 100.0 * k / len(st), 100.0 * k / len(hits), note))
        print("%-36s size %6d  AntiFam hits %4d (%.2f%% of set)" % (nm, len(st), k, 100.0 * k / len(st)))
with open(OUT2, "w") as o:
    o.write("rep\tantifam_families\t" + "\t".join(nm for nm, _, _ in S[1:]) + "\n")
    for r in sorted(hits): o.write("%s\t%s\t%s\n" % (r, ",".join(sorted(hitfam[r])), "\t".join("1" if r in st else "0" for _, st, _ in S[1:])))
print("wrote", OUT1, OUT2, "and", ", ".join(LISTS))
