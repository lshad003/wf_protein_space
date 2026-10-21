# Step 34e summary: old whole-catalog n88 values beside the new strict_final and strict_final EU values.
# New: results/step34e_{strict_final,strict_final_EU}_{permanova,limma}.tsv (step34e_hostfactor.R).
# Old (whole n88 catalog, 2,067,011 supported reps): logs/step22b_n88.log, logs/step22c_n88.log, logs/step22e_n88.log,
#   logs/step17g_n88.log, logs/step17h_n88.log, logs/step17i_n88.log (values copied below from those logs).
# Output: results/step34e_summary.tsv; refuses to overwrite.
# Optional argument (Step 34s): <dir>; reads the new tables from <dir> and writes <dir>/step34e_summary.tsv. The old
#   columns stay the whole-catalog values below. Default (no argument) unchanged.
import os, sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
D = sys.argv[1].rstrip("/") + "/" if len(sys.argv) > 1 else R
OUT = D + "step34e_summary.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
A, P = "abundance_bray_hellinger", "presence_jaccard_ge5"
OLDP = {("WF22", A, "whole_animal_perm_balanced_42", "egg_mass"): ("0.32728", "0.001", "step22b_n88"),
        ("WF22", A, "whole_animal_perm_balanced_42", "treatment"): ("0.05121", "0.928", "step22b_n88"),
        ("WF22", A, "within_animal_perm_44", "month"): ("0.12698", "0.001", "step22b_n88"),
        ("WF22", A, "animal_centroids_15", "egg_mass"): ("0.51267", "0.001", "step22b_n88"),
        ("WF22", A, "animal_centroids_15", "treatment"): ("0.06948", "0.657", "step22b_n88"),
        ("WF22", A, "betadisper_centroids", "treatment"): ("F 9.6328", "0.003197", "step22b_n88"),
        ("WF22", A, "betadisper_centroids", "egg_mass"): ("F 1.7605", "0.2136", "step22b_n88"),
        ("WF22", P, "whole_animal_perm_balanced_42", "egg_mass"): ("0.15651", "0.001", "step22c_n88"),
        ("WF22", P, "whole_animal_perm_balanced_42", "treatment"): ("0.06205", "0.665", "step22c_n88"),
        ("WF22", P, "animal_centroids_15", "egg_mass"): ("0.36301", "0.001", "step22c_n88"),
        ("WF22", P, "animal_centroids_15", "treatment"): ("0.10570", "0.465", "step22c_n88"),
        ("WF22", P, "richness_kruskal_44", "treatment"): ("chi2 2.0014", "0.3676", "step22c_n88"),
        ("WF24_n34", A, "free_perm", "Model"): ("0.17251", "0.679", "step22e_n88"),
        ("WF24_n34", P, "free_perm", "Model"): ("0.16957", "0.77", "step22e_n88")}
OLDL = {("WF22", "B_naive", "treatmentSTP1710.7"): 5005, ("WF22", "B_naive", "treatmentSTP1717.1"): 7263,
        ("WF22", "B_naive", "egg_mass2"): 84645, ("WF22", "B_naive", "egg_mass3"): 69732,
        ("WF22", "B_naive", "month2"): 68938, ("WF22", "B_naive", "month3"): 117641,
        ("WF22", "A_blocked", "within_animal_correlation"): 0.2864,
        ("WF22", "A_blocked", "treatmentSTP1710.7"): 3, ("WF22", "A_blocked", "treatmentSTP1717.1"): 155,
        ("WF22", "A_blocked", "egg_mass2"): 57120, ("WF22", "A_blocked", "egg_mass3"): 51808,
        ("WF22", "A_blocked", "month2"): 99067, ("WF22", "A_blocked", "month3"): 141580,
        ("WF22", "C_EM3", "treatmentSTP1710.7"): 0, ("WF22", "C_EM3", "month2"): 80470, ("WF22", "C_EM3", "month3"): 110406,
        ("WF23", "limma", "treatmentUHM520.7734"): 0}
for t, n in (("1", 0), ("2", 0), ("3", 3), ("5", 0), ("6", 0), ("7", 0)): OLDL[("WF24_n34", "limma", "treatment" + t)] = n
OLDG = {"WF22": 367923, "WF23": 565822, "WF24_n34": 493896}
def rd(p):
    rows = [l.rstrip("\n").split("\t") for l in open(p)]; return [dict(zip(rows[0], r)) for r in rows[1:]]
new = {s: (rd(D + "step34e_%s_permanova.tsv" % s), rd(D + "step34e_%s_limma.tsv" % s)) for s in ("strict_final", "strict_final_EU")}
out = [("kind", "cohort", "metric_or_model", "test", "term_or_coef", "old_whole_n88", "old_p", "old_source",
        "strict_final", "strict_final_p", "strict_final_EU", "strict_final_EU_p", "note")]
keys = [(r["cohort"], r["metric"], r["test"], r["term"]) for r in new["strict_final"][0]]
for k in keys:
    a = [r for r in new["strict_final"][0] if (r["cohort"], r["metric"], r["test"], r["term"]) == k][0]
    e = [r for r in new["strict_final_EU"][0] if (r["cohort"], r["metric"], r["test"], r["term"]) == k][0]
    val = lambda r: ("F " + r["F"]) if r["test"] in ("betadisper_centroids",) else (("chi2 " + r["F"]) if r["test"].startswith("richness") else r["R2"])
    o = OLDP.get(k, ("NA", "NA", "not run"))
    out.append(("permanova",) + k + (o[0], o[1], o[2], val(a), a["p"], val(e), e["p"], a["note"]))
lk = [(r["cohort"], r["model"], r["coef"]) for r in new["strict_final"][1]]
for k in lk:
    a = [r for r in new["strict_final"][1] if (r["cohort"], r["model"], r["coef"]) == k][0]
    e = [r for r in new["strict_final_EU"][1] if (r["cohort"], r["model"], r["coef"]) == k][0]
    o = OLDL.get(k, "NA"); g = OLDG.get(k[0], "NA")
    note = "genes tested: old %s, strict_final %s, EU %s" % (g, a["genes_tested"], e["genes_tested"])
    out.append(("limma",) + (k[0], k[1], "padj<0.05 (correlation for within_animal_correlation)", k[2], str(o), "", "step17g/h/i_n88" if o != "NA" else "not run",
               a["padj_lt_0.05"], "", e["padj_lt_0.05"], "", note))
with open(OUT, "w") as fh:
    for r in out: fh.write("\t".join(map(str, r)) + "\n")
for r in out: print("\t".join(map(str, r)))
print("wrote", OUT)
