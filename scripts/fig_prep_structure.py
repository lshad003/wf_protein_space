# Prep for Figure 4. 1) Best ProstT5/Foldseek E-value (min over PDB and CATH50) per sequence: strict_final recurrent EU
# family reps (queries), fixed length-matched K controls, reversed decoys (Steps 34d/34d2); no hit = NA.
# 2) Prodigal coordinates of the genes on contig UHM636.41030__k141_332084 (Cse2/CasB-like), with the Pfam labels of
# the +-3 neighbours of member _5 from step29d. Outputs: results/fig_prep_bestE.tsv, results/fig_prep_cse2_locus.tsv
import csv, glob, subprocess
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
def fasta(p):
    L, cur = {}, None
    for l in open(p):
        if l.startswith(">"): cur = l[1:].split()[0]; L[cur] = 0
        else: L[cur] += len(l.strip().rstrip("*"))
    return L
SETS = {"queries": fasta(R + "step34d_queries.faa"), "K": fasta(R + "step34d_known_lenmatched_fixed.faa"), "decoys": fasta(R + "step34d_decoys_rev.faa")}
print("sequences:", {k: len(v) for k, v in SETS.items()})
m8 = sorted(glob.glob(R + "step34d_prostt5/*_vs_*.tsv") + glob.glob(R + "step34d2_prostt5/*_vs_*.tsv"))
print("m8 files:", [f.split("/")[-1] for f in m8]); best = {}
for f in m8:
    for l in open(f):
        p = l.split("\t")
        try: e = float(p[2])
        except (IndexError, ValueError): continue
        if p[0] not in best or e < best[p[0]]: best[p[0]] = e
with open(R + "fig_prep_bestE.tsv", "w") as o:
    o.write("set\tid\tlength\tbestE\n")
    for s, L in SETS.items():
        nh = 0
        for i, n in L.items():
            e = best.get(i); nh += (e is not None and e <= 1e-3)
            o.write("%s\t%s\t%d\t%s\n" % (s, i, n, "NA" if e is None else "%.6g" % e))
        print("%s: %d sequences, %d with best E <= 1e-3 (expected queries 186, K 2264, decoys 15)" % (s, len(L), nh))
C = "UHM636.41030__k141_332084"
hd = subprocess.run("zcat %s/catalog/db/LsPS_AA_all.fasta.gz | grep '^>%s_'" % (W, C), shell=True, capture_output=True, text=True).stdout.splitlines()
lab = {}
for r in csv.DictReader(open(R + "step29d_fold4_neighbours.tsv"), delimiter="\t"):
    if r["member"] == C + "_5" and r["pfam"] not in ("", "NA", "-"): lab.setdefault(r["neighbour"], set()).add(r["pfam"])
with open(R + "fig_prep_cse2_locus.tsv", "w") as o:
    o.write("gene\tstart\tend\tstrand\tlabel\n")
    for h in hd:
        f = h[1:].split(" # "); gid = f[0].split()[0]
        lb = "this study" if gid == C + "_5" else ",".join(sorted(lab.get(gid, ()))) or "-"
        o.write("%s\t%s\t%s\t%s\t%s\n" % (gid, f[1], f[2], f[3], lb)); print(gid, f[1], f[2], f[3], lb)
