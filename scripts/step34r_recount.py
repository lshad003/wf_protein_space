# Step 34r: read-based numbers recomputed on the identity-filtered matrices, side by side with the old matrix
#   (results/count_matrix_primary_n88.tsv.gz). Old matrix is run through the same code as a reproduction check.
# Matrices: old, col2, c_id95, c_id90, c_id95_q1 (results/step34r/matrix_*_n88.tsv.gz, scripts/step34r_matrices.py).
# 1. Prevalence, rules of scripts/step34b2_final_stats.py: strict_final reps (results/step34a_rep_origin_v2_final.tsv
#    strict_final = 1), class from results/fourway_classes_n88.tsv; EU 41,114, K all strict_final K (1,118,691).
#    Mean number of the 88 metagenomes with >= 1 and >= 5 reads; % of reps with >= 1 (>= 5) reads in a sample of each of
#    the three cohorts (metadata/samples_88_treatment.tsv).
# 2. WF22 egg-mass D at >= 5 reads, rules of scripts/step34c3b_checks.py: WF22 months 2-3 (15 animals), presence = >= 5
#    reads in any retained sample; 41,114 EU and 41,114 random strict_final K reps (seed 343); restricted = >= 3 and
#    >= 50% of the egg mass's animals and absent from the other egg masses; D = EU - K summed over egg masses; null from
#    1,000 permutations of egg-mass labels (numpy default_rng(34)), EU and K together; p = (1 + #null >= obs) / 1,001.
# 3. WF24 egg mass 4 tank-restricted EU vs K at >= 5 reads, rules of scripts/step34m_tank_restriction.py: tanks 4.1 / 4.2
#    (metadata/wf23_wf24_sequenced_metadata.tsv), UHM586.41010 and UHM590.41012 excluded (metadata/wf24_excluded.tsv),
#    one sample per animal, same units and restriction rule, both tanks summed, fresh default_rng(34), 1,000 permutations.
# Old values that must reproduce on the old matrix: EU 37.05 / 22.72, 72.13% / 35.17%; D 454; WF24 EU 723, K 284.
# Outputs: results/step34r/recount_n88.tsv (long), results/step34r/recount_side_by_side_n88.tsv; refuse to overwrite.
import os, random
import numpy as np, pandas as pd
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; O = R + "step34r/"
OUT1, OUT2 = O + "recount_n88.tsv", O + "recount_side_by_side_n88.tsv"
for p in (OUT1, OUT2):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
MATS = [("old", R + "count_matrix_primary_n88.tsv.gz")] + [(n, O + "matrix_%s_n88.tsv.gz" % n) for n in ("col2", "c_id95", "c_id90", "c_id95_q1")]
NPERM = 1000
meta = pd.read_csv(W + "/metadata/samples_88_treatment.tsv", sep="\t", dtype=str).set_index("stem")
orig = pd.read_csv(R + "step34a_rep_origin_v2_final.tsv", sep="\t", dtype=str, usecols=["rep", "class4", "strict_final"])
fw = pd.read_csv(R + "fourway_classes_n88.tsv", sep="\t", dtype=str, usecols=["rep", "class4"]).set_index("rep")["class4"]
fin = orig[orig.strict_final == "1"]
EU = sorted(fin.rep[fin.class4 == "EU"]); Kc = sorted(fin.rep[fin.class4 == "K"])
if len(EU) != 41114: raise SystemExit("ABORT: expected 41,114 EU")
if (fw[EU] != "EU").any() or (fw[Kc] != "K").any(): raise SystemExit("ABORT: class4 differs between tables")
random.seed(343); Kb = sorted(random.sample(sorted(Kc), len(EU)))
# WF22 design (34c3b)
d22 = pd.read_csv(W + "/metadata/wf22_design.tsv", sep="\t", dtype=str).set_index("stem")
keep22 = [s for s in meta.index if meta.cohort[s] == "WF22" and int(d22.month[s]) in (2, 3)]
an22 = sorted({s.split(".")[0] for s in keep22}); em22 = {s.split(".")[0]: meta.egg_mass_cohort[s] for s in keep22}
lab22 = [em22[a] for a in an22]
# WF24 EM4 tanks (34m)
excl = set(pd.read_csv(W + "/metadata/wf24_excluded.tsv", sep="\t", dtype=str).iloc[:, 0])
if excl != {"UHM586.41010", "UHM590.41012"}: raise SystemExit("ABORT: exclusions %s" % excl)
tm = pd.read_csv(W + "/metadata/wf23_wf24_sequenced_metadata.tsv", sep="\t", dtype=str).set_index("stem")
st24 = [s for s in meta.index if s not in excl and meta.cohort[s] != "WF22" and tm.wf_year[s] + "_EM" + tm.egg_mass[s] == "WF24_EM4"]
an24 = sorted(s.split(".")[0] for s in st24); s24 = {s.split(".")[0]: s for s in st24}
if len(an24) != len(st24): raise SystemExit("ABORT: WF24 EM4 animal with two samples")
lab24 = [tm.tank[s24[a]] for a in an24]
print("WF22 animals %d samples %d; WF24 EM4 animals %d tanks %s" % (len(an22), len(keep22), len(an24),
      {t: lab24.count(t) for t in sorted(set(lab24))}), flush=True)
def restricted(M, lab):
    tot = M.sum(1); out = {}
    for e in sorted(set(lab)):
        g = np.array([x == e for x in lab]); n = M[:, g].sum(1)
        out[e] = (n >= 3) & (n >= 0.5 * g.sum()) & (tot - n == 0)
    return out
def paired(PE, PK, lab0):
    G = sorted(set(lab0)); rng = np.random.default_rng(34); perms = [list(rng.permutation(lab0)) for _ in range(NPERM)]
    def score(lab):
        cE, cK = restricted(PE, lab), restricted(PK, lab)
        return sum(int(cE[g].sum()) for g in G), sum(int(cK[g].sum()) for g in G)
    oE, oK = score(lab0); N = np.array([score(l) for l in perms]); nd = N[:, 0] - N[:, 1]
    p = lambda o, n: (1 + int((n >= o).sum())) / (NPERM + 1)
    return {"EU": (oE, N[:, 0].mean(), N[:, 0].std(), p(oE, N[:, 0])), "K": (oK, N[:, 1].mean(), N[:, 1].std(), p(oK, N[:, 1])),
            "D": (oE - oK, nd.mean(), nd.std(), p(oE - oK, nd))}
rows = []
def add(m, metric, cls, val, null_mean="", null_sd="", p=""):
    rows.append((m, metric, cls, val, null_mean, null_sd, p)); print("  %-10s %-34s %-3s %s %s %s %s" % (m, metric, cls, val, null_mean, null_sd, p), flush=True)
for m, f in MATS:
    print("reading", m, f, flush=True)
    df = pd.read_csv(f, sep="\t", index_col=0, dtype={0: str})
    if len(df) != 2067011 or list(df.columns) != list(pd.read_csv(MATS[0][1], sep="\t", nrows=0, index_col=0).columns):
        raise SystemExit("ABORT: matrix layout " + m)
    coh = np.array([meta.cohort[s] for s in df.columns])
    for cls, L in (("EU", EU), ("K", Kc)):
        X = df.loc[L].values
        for t in (1, 5):
            B = X >= t
            add(m, "mean_prevalence_ge%d" % t, cls, "%.2f" % B.sum(1).mean())
            c3 = np.all([B[:, coh == c].any(1) for c in ("WF22", "WF23", "WF24")], axis=0)
            add(m, "pct_3cohorts_ge%d" % t, cls, "%.2f" % (100.0 * c3.mean()))
    def pres(L, stems_by_animal):
        X = df.loc[L]
        return np.stack([(X[ss].values >= 5).any(1) for ss in stems_by_animal], 1)
    sb22 = [[s for s in keep22 if s.split(".")[0] == a] for a in an22]
    r = paired(pres(EU, sb22), pres(Kb, sb22), lab22)
    for k in ("EU", "K", "D"):
        o, nm, sd, p = r[k]; add(m, "WF22_eggmass_restricted_ge5", k, o, "%.1f" % nm, "%.1f" % sd, "%.4f" % p)
    sb24 = [[s24[a]] for a in an24]
    r = paired(pres(EU, sb24), pres(Kb, sb24), lab24)
    for k in ("EU", "K", "D"):
        o, nm, sd, p = r[k]; add(m, "WF24_EM4_tank_restricted_ge5", k, o, "%.1f" % nm, "%.1f" % sd, "%.4f" % p)
    del df
with open(OUT1, "w") as o:
    o.write("matrix\tmetric\tclass\tvalue\tnull_mean\tnull_sd\tp_one_sided\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
t = pd.DataFrame(rows, columns=["matrix", "metric", "class", "value", "null_mean", "null_sd", "p"])
sbs = t.pivot_table(index=["metric", "class"], columns="matrix", values="value", aggfunc="first")[[m for m, _ in MATS]]
sbs.to_csv(OUT2, sep="\t")
print(sbs.to_string())
old = {("mean_prevalence_ge1", "EU"): "37.05", ("mean_prevalence_ge5", "EU"): "22.72", ("pct_3cohorts_ge1", "EU"): "72.13",
       ("pct_3cohorts_ge5", "EU"): "35.17", ("WF22_eggmass_restricted_ge5", "D"): "454",
       ("WF24_EM4_tank_restricted_ge5", "EU"): "723", ("WF24_EM4_tank_restricted_ge5", "K"): "284"}
bad = [(k, v, sbs.loc[k, "old"]) for k, v in old.items() if str(sbs.loc[k, "old"]) != v]
print("old values reproduced on the old matrix:", "yes" if not bad else "NO %s" % bad)
print("wrote", OUT1, OUT2)
if bad: raise SystemExit(1)
