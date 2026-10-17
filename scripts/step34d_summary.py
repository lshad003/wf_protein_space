# Step 34d summary: ProstT5/Foldseek hits of the strict_final recurrent EU queries, length-matched K and reversed decoys.
# Hit = best E <= 1e-3 against PDB or CATH50 (min over both; column 3), as step28h. Lengths without '*'; bins
#   <100, 100-199, 200-399, >=400 aa (as step28g2/28h). Decoy false-positive rate = decoy hits / decoys, overall and per
#   bin; expected false hits among query hits = FPR x queries (per bin, summed).
# Old HIT: results/prostt5_main_hit_nostar.tsv status kept or gained (214), sequences from
#   results/prostt5_main_eu15580_nostar.faa; matched to the new queries by exact sequence.
# Cse2 family rep UHM45.35758__k141_100597_236.
# Outputs (refuse to overwrite): results/step34d_by_length.tsv, results/step34d_query_hits.tsv
import os
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34d_prostt5/"
O1, O2 = R + "step34d_by_length.tsv", R + "step34d_query_hits.tsv"
for p in (O1, O2):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
EMAX = 1e-3; CSE2 = "UHM45.35758__k141_100597_236"
BINS = (("<100", 0, 100), ("100-199", 100, 200), ("200-399", 200, 400), (">=400", 400, 10 ** 9), ("all", 0, 10 ** 9))
def read(faa):
    s, h = {}, None
    for l in open(faa):
        if l[0] == ">": h = l[1:].split()[0]; s[h] = []
        else: s[h].append(l.strip())
    return {k: "".join(v) for k, v in s.items()}
def best(pre, ids):
    b = {}; tg = {}
    for db in ("PDB", "CATH50"):
        for l in open(D + "%s_vs_%s.tsv" % (pre, db)):
            f = l.split("\t")
            if f[0] not in ids: raise SystemExit("ABORT: unknown query in %s_vs_%s.tsv" % (pre, db))
            e = float(f[2])
            if e < b.get(f[0], 9e9): b[f[0]] = e; tg[f[0]] = "%s:%s" % (db, f[1])
    return b, tg
S = {}
for lab, faa, pre in (("queries", "step34d_queries.faa", "queries"), ("K_lenmatched", "step34d_known_lenmatched.faa", "known"),
                      ("decoys_rev", "step34d_decoys_rev.faa", "decoys")):
    sq = read(R + faa); b, tg = best(pre, sq)
    S[lab] = (sq, {k for k, e in b.items() if e <= EMAX}, b, tg)
    print("%s: %d sequences, with any row %d, hit (E <= 1e-3) %d" % (lab, len(sq), len(b), len(S[lab][1])))
rows = [("set", "len_bin_aa", "n", "n_hit", "pct_hit")]
fpr = {}
print("\nhit rate by length bin (best E <= 1e-3, PDB or CATH50)")
print("%-9s %22s %22s %22s %18s" % ("bin", "queries", "K_lenmatched", "decoys_rev", "expected false hits"))
exp_tot = 0.0
for bn, lo, hi in BINS:
    out = []
    for lab in ("queries", "K_lenmatched", "decoys_rev"):
        sq, H = S[lab][0], S[lab][1]
        ids = [k for k, v in sq.items() if lo <= len(v) < hi]; n = len(ids); h = sum(1 for k in ids if k in H)
        rows.append((lab, bn, n, h, "%.2f" % (100.0 * h / n) if n else "NA")); out.append((n, h))
        if lab == "decoys_rev": fpr[bn] = h / n if n else 0.0
    ex = fpr[bn] * out[0][0]
    if bn != "all": exp_tot += ex
    print("%-9s %6d/%-6d (%5.2f%%) %6d/%-6d (%5.2f%%) %6d/%-6d (%5.2f%%) %18.1f" % (bn, out[0][1], out[0][0], 100.0 * out[0][1] / max(1, out[0][0]),
          out[1][1], out[1][0], 100.0 * out[1][1] / max(1, out[1][0]), out[2][1], out[2][0], 100.0 * out[2][1] / max(1, out[2][0]), ex))
nq = len(S["queries"][0]); hq = len(S["queries"][1])
print("decoy FPR overall %.4f%%; expected false hits among %d query hits: %.1f (overall FPR x queries) / %.1f (summed by bin)"
      % (100.0 * fpr["all"], hq, fpr["all"] * nq, exp_tot))
with open(O1, "w") as o:
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
qs, qh, qb, qt = S["queries"]
old = {}
for i, l in enumerate(open(R + "prostt5_main_hit_nostar.tsv")):
    if i: f = l.rstrip("\n").split("\t"); old[f[0]] = f[4]
o214 = {k for k, s in old.items() if s in ("kept", "gained")}
if len(o214) != 214: raise SystemExit("ABORT: old HIT is %d, expected 214" % len(o214))
oseq = read(R + "prostt5_main_eu15580_nostar.faa")
bys = {}
for k, v in qs.items(): bys.setdefault(v, []).append(k)
m = {k: bys.get(oseq[k], []) for k in o214}
inq = [k for k in o214 if m[k]]
print("\nold 214 no-star HIT families: %d have their sequence among the new queries; of those, HIT again %d" % (len(inq), sum(1 for k in inq if any(x in qh for x in m[k]))))
lab = {}
for k in o214:
    for x in m[k]: lab[x] = k
with open(O2, "w") as o:
    o.write("query\tlength_aa\tbest_E\tbest_target\thit\told214_rep_same_sequence\n")
    for k in sorted(qs):
        e = qb.get(k)
        o.write("%s\t%d\t%s\t%s\t%d\t%s\n" % (k, len(qs[k]), "%.3g" % e if e is not None else "NA", qt.get(k, "NA"), int(k in qh), lab.get(k, "-")))
print("Cse2 rep %s: in queries %s | best E %s (%s) | hit %s (expected no hit)" % (CSE2, CSE2 in qs, qb.get(CSE2, "none"), qt.get(CSE2, "-"), CSE2 in qh))
print("wrote", O1, "and", O2)
