# Prep for Figure 3. Level 2: per-rep within-year recurrence flags (assembly presence, one sample per animal, WF22
# month MONTH, UHM590 excluded; same EU set and K sample as step34c4). Level 3: presence across WF22 animals (months 2
# and 3, as step34c3) of the egg-mass-restricted EU reps. Outputs: results/fig_prep_year_flags.tsv,
# results/fig_prep_wf22_em_presence.tsv
import csv, random
import numpy as np
MONTH = "2"
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
eu, kpool = [], []
with open(R + "step34a_rep_origin_v2_final.tsv") as f:
    h = f.readline().rstrip("\n").split("\t"); ir, ic, isf = h.index("rep"), h.index("class4"), h.index("strict_final")
    for l in f:
        p = l.rstrip("\n").split("\t")
        if p[isf] in ("1", "True", "true", "yes"):
            if p[ic] == "EU": eu.append(p[ir])
            elif p[ic] == "K": kpool.append(p[ir])
random.seed(3443); kk = random.sample(sorted(kpool), len(eu)); SETS = {"EU": eu, "K": kk}; target = set(eu) | set(kk)
meta = {r["stem"]: r for r in csv.DictReader(open(W + "/metadata/samples_88_treatment.tsv"), delimiter="\t") if not r["stem"].startswith("UHM590.")}
des = {r["stem"]: r for r in csv.DictReader(open(W + "/metadata/wf22_design.tsv"), delimiter="\t")}
pres = {}
with open(R + "clusters_95_n88.tsv") as f:
    for l in f:
        p = l.rstrip("\n").split("\t")
        if p[0] in target: pres.setdefault(p[0], set()).add(p[1].split("__")[0])
stems = [s for s in meta if meta[s]["cohort"] != "WF22" or des.get(s, {}).get("month") == MONTH]; ss = set(stems)
acoh = {s.split(".")[0]: meta[s]["cohort"] for s in stems}; Y = ["WF22", "WF23", "WF24"]
thr = {y: max(3, int(np.ceil(0.5 * sum(1 for a in acoh if acoh[a] == y)))) for y in Y}
out = open(R + "fig_prep_year_flags.tsv", "w"); out.write("set\trep\twr_WF22\twr_WF23\twr_WF24\n")
for sn, reps in SETS.items():
    tally = dict.fromkeys(Y, 0); core = 0
    for rp in reps:
        an = {st.split(".")[0] for st in pres.get(rp, ()) if st in ss}
        wr = [int(sum(1 for a in an if acoh[a] == y) >= thr[y]) for y in Y]; out.write("%s\t%s\t%d\t%d\t%d\n" % (sn, rp, *wr))
        for y, v in zip(Y, wr): tally[y] += v
        core += all(wr)
    print("level 2 %s: within-year recurrent %s, core %d" % (sn, tally, core))
out.close()
w22 = [s for s in des if des[s]["month"] in ("2", "3")]; s22 = set(w22)
an22 = sorted({des[s]["animal"] for s in w22}); em = {des[s]["animal"]: des[s]["egg_mass"] for s in w22}
trt = {des[s]["animal"]: des[s]["treatment"] for s in w22}; EM = sorted(set(em.values()))
nem = {e: sum(1 for a in an22 if em[a] == e) for e in EM}; cnt = dict.fromkeys(EM, 0)
o = open(R + "fig_prep_wf22_em_presence.tsv", "w")
o.write("rep\trestricted_em\t" + "\t".join("%s|%s|%s" % (a, em[a], trt[a]) for a in an22) + "\n")
for rp in eu:
    an = {st.split(".")[0] for st in pres.get(rp, ()) if st in s22}; c = {e: sum(1 for a in an if em[a] == e) for e in EM}
    for e in EM:
        if c[e] >= max(3, int(np.ceil(0.5 * nem[e]))) and sum(c[x] for x in EM if x != e) == 0:
            cnt[e] += 1; o.write("%s\tEM%s\t%s\n" % (rp, e, "\t".join("1" if a in an else "0" for a in an22)))
o.close(); print("level 3 EU restricted per egg mass:", cnt, "(Step 34c3: 48 / 63 / 837)")
