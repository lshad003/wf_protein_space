# Step 34g Part C: which Cse2 variant Figure 4c (ESMFold model) and Figure 4d (rank-1 locus) show.
# Inputs: results/esmfold_pdb/UHM636.41030__k141_332084_5.pdb (sequence from CA records), catalog proteins via
#   results/step34f_cse2/prodigal.faa (identical to catalog/db/LsPS_AA_all.fasta.gz for all 10,476 genes, Step 34f),
#   results/step34f_cse2/members.json (family members), results/fig_prep_cse2_locus_best.tsv (rank-1 locus genes),
#   results/step34f_cse2/cse2_nocut.tblout (CRISPR_Cse2 PF09485.16 scores without thresholds, all contig genes),
#   results/step34f_cse2/hmmsearch.tblout (--cut_ga). Variant A = the full-length sequence of the rep
#   UHM45.35758__k141_100597_236; variant B = the most frequent other full-length (161 aa) member sequence.
# Output: results/step34g_cse2_model_check.tsv (refuses to overwrite)
import os, json, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; F = R + "step34f_cse2/"
OUT = R + "step34g_cse2_model_check.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
AA = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I", LEU="L", LYS="K",
          MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
pdb = "".join(AA[l[17:20]] for l in open(R + "esmfold_pdb/UHM636.41030__k141_332084_5.pdb") if l.startswith("ATOM") and l[12:16].strip() == "CA")
seq = {}; h = None
for l in open(F + "prodigal.faa"):
    if l[0] == ">": h = l[1:].split()[0]; seq[h] = ""
    else: seq[h] += l.strip()
seq = {k: v.rstrip("*") for k, v in seq.items()}
mem = json.load(open(F + "members.json"))["mem"]
REP = "UHM45.35758__k141_100597_236"; MOD = "UHM636.41030__k141_332084_5"
best = [l.split("\t") for l in open(R + "fig_prep_cse2_locus_best.tsv")]
loc = [r for r in best if r[0] == "gene" and "cse2" in r[7].split(",")]
LOC = loc[0][1]
A = seq[REP]
full = collections.Counter(seq[m] for m in mem if len(seq[m]) == len(A))
B = [s for s, n in full.most_common() if s != A][0]
nA = sum(1 for m in mem if seq[m] == A); nB = sum(1 for m in mem if seq[m] == B)
def cmp(x, y):
    if len(x) != len(y): return "NA", "length %d vs %d" % (len(x), len(y))
    d = [(i + 1, a, b) for i, (a, b) in enumerate(zip(x, y)) if a != b]
    return "%.2f" % (100.0 * (len(x) - len(d)) / len(x)), ";".join("%s%d%s" % (a, i, b) for i, a, b in d) or "none"
nocut = {}
for l in open(F + "cse2_nocut.tblout"):
    if l[0] != "#": p = l.split(); nocut[p[0]] = (float(p[5]), p[4])
ga = {l.split()[0] for l in open(F + "hmmsearch.tblout") if l[0] != "#" and l.split()[2] == "CRISPR_Cse2"}
GA = 24.6
rows = [("check", "value")]
rows.append(("model sequence length (CA records)", len(pdb)))
rows.append(("model identical to catalog protein " + MOD, pdb == seq[MOD]))
rows.append(("variant A = rep " + REP + ", members", nA))
rows.append(("variant B = most frequent other 161 aa sequence, members", nB))
rows.append(("model identical to variant A", pdb == A)); rows.append(("model identical to variant B", pdb == B))
for lab, s in (("rep / variant A", A), ("variant B", B), ("rank-1 locus cse2 " + LOC, seq[LOC])):
    pid, dif = cmp(pdb, s); rows.append(("model vs %s: percent identity" % lab, pid)); rows.append(("model vs %s: differences (model pos other)" % lab, dif))
pid, dif = cmp(A, B); rows.append(("variant A vs variant B: percent identity", pid)); rows.append(("variant A vs variant B: differences", dif))
rows.append(("rank-1 locus cse2 identical to variant A", seq[LOC] == A)); rows.append(("rank-1 locus cse2 identical to variant B", seq[LOC] == B))
for lab, g in (("model " + MOD, MOD), ("rep " + REP, REP), ("rank-1 locus cse2 " + LOC, LOC)):
    sc = nocut.get(g)
    rows.append(("%s: CRISPR_Cse2 bit score without --cut_ga (E)" % lab, "%.1f (%s)" % sc if sc else "no row at T 0"))
    rows.append(("%s: hit with --cut_ga (GA %.1f)" % (lab, GA), g in ga))
    if sc: rows.append(("%s: margin to GA (bits)" % lab, "%.1f" % (sc[0] - GA)))
L = {r[1]: r for r in best if r[0] == "gene"}
cas = collections.Counter(n for r in L.values() for n in r[7].split(",") if n != "-")
rows.append(("rank-1 locus contig", LOC.rsplit("_", 1)[0]))
for n in ("cas3", "cse1", "cse2", "cas7", "cas5", "cas6e", "cas1", "cas2"): rows.append(("rank-1 locus carries " + n, cas[n] > 0))
rows.append(("rank-1 locus CRISPR arrays (MinCED 0.4.2)", sum(1 for r in best if r[0] == "array")))
with open(OUT, "w") as o:
    for r in rows: o.write("%s\t%s\n" % r)
for r in rows: print("%s\t%s" % r)
