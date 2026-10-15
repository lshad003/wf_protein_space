# Step 28h compare: original Step 28 searches (query ends in '*') against the _nostar reruns (step28h_nostar.sh).
# Sets (original FASTA and hits -> nostar FASTA and hits):
#  EU15580          results/prostt5_main/eu15580_vs_{db}.tsv           -> results/prostt5_main_nostar/eu15580_vs_{db}.tsv
#  EU740            results/prostt5_pilot/eu740_vs_{db}.tsv (28e2)     -> results/prostt5_pilot_nostar/eu740_vs_{db}.tsv
#  K740_lenmatched  results/prostt5_calib/known740_vs_{db}.tsv         -> results/prostt5_calib_nostar/known740_vs_{db}.tsv
#  EU740_reversed   results/prostt5_calib/decoy740_vs_{db}.tsv         -> results/prostt5_calib_nostar/decoy740_vs_{db}.tsv
# Best E per query from column 3; db = PDB, CATH50, either (min over both). Lengths = residues without '*'
#   (as step28g2), so bins are the same in both runs. Denominator: sequences in each FASTA.
# Outputs (refuse to overwrite):
#  results/prostt5_thresholds_nostar.tsv  queries with a hit at E <= 1e-1..1e-10 per set x db x run (as 28f/28g summaries)
#  results/prostt5_by_length_nostar.tsv   hit (either, E <= 1e-3) by length bin per set, original vs nostar (as 28g2)
#  results/prostt5_main_hit_nostar.tsv    EU15580 reps that are HIT (either, E <= 1e-3) in either run: best E in each
#                                          run and status kept / lost / gained. Aborts unless the original HIT is 218.
import os
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
O1, O2, O3 = R + "prostt5_thresholds_nostar.tsv", R + "prostt5_by_length_nostar.tsv", R + "prostt5_main_hit_nostar.tsv"
for p in (O1, O2, O3):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
EMAX = 1e-3; TH = (1e-1, 1e-2, 1e-3, 1e-5, 1e-10)
SETS = (("EU15580", "prostt5_main_eu15580", "prostt5_main", "eu15580", 15580),
        ("EU740", "prostt5_pilot_eu740", "prostt5_pilot", "eu740", 740),
        ("K740_lenmatched", "prostt5_calib_known_lenmatched", "prostt5_calib", "known740", 740),
        ("EU740_reversed", "prostt5_decoy_eu740rev", "prostt5_calib", "decoy740", 740))
BINS = (("<100", 0, 100), ("100-199", 100, 200), ("200-399", 200, 400), (">=400", 400, 10 ** 9), ("all", 0, 10 ** 9))
def lengths(faa):
    L, h = {}, None
    for l in open(faa):
        if l[0] == ">": h = l[1:].split()[0]; L[h] = 0
        else: L[h] += len(l.strip().replace("*", ""))
    return L
def best(tsv, q, qlen):
    b = {}
    for l in open(tsv):
        f = l.split("\t")
        if f[0] not in q: raise SystemExit("ABORT: hit query not in query set: " + tsv)
        e = float(f[2]); b[f[0]] = min(e, b.get(f[0], e)); qlen[f[0]] = int(f[6])
    return b
t1 = [("set", "db", "run", "n_queries") + tuple("n_E<=%g" % t for t in TH) + tuple("pct_E<=%g" % t for t in TH)]
t2 = [("set", "len_bin_aa", "n_queries", "n_hit_orig", "pct_hit_orig", "n_hit_nostar", "pct_hit_nostar", "n_kept", "n_lost", "n_gained")]
B = {}
for lab, faa, d, pre, n in SETS:
    L = lengths(R + faa + ".faa"); L2 = lengths(R + faa + "_nostar.faa")
    if len(L) != n or L != L2: raise SystemExit("ABORT: %s: %d sequences or lengths differ from the nostar FASTA" % (lab, len(L)))
    ql = {"orig": {}, "nostar": {}}
    for run, dd in (("orig", d), ("nostar", d + "_nostar")):
        b = {db: best(R + "%s/%s_vs_%s.tsv" % (dd, pre, db), L, ql[run]) for db in ("PDB", "CATH50")}
        b["either"] = {x: min(b["PDB"].get(x, 9e9), b["CATH50"].get(x, 9e9)) for x in set(b["PDB"]) | set(b["CATH50"])}
        B[(lab, run)] = b
        for db in ("PDB", "CATH50", "either"):
            k = [sum(1 for e in b[db].values() if e <= t) for t in TH]
            t1.append((lab, db, run, n) + tuple(k) + tuple("%.2f" % (100.0 * x / n) for x in k))
    H = {run: {x for x, e in B[(lab, run)]["either"].items() if e <= EMAX} for run in ("orig", "nostar")}
    for bn, lo, hi in BINS:
        ids = {x for x in L if lo <= L[x] < hi}; a, c = H["orig"] & ids, H["nostar"] & ids
        pc = lambda s: "%.2f" % (100.0 * len(s) / len(ids)) if ids else "NA"
        t2.append((lab, bn, len(ids), len(a), pc(a), len(c), pc(c), len(a & c), len(a - c), len(c - a)))
    both = set(ql["orig"]) & set(ql["nostar"]); dq = {}
    for x in both: k = ql["orig"][x] - ql["nostar"][x]; dq[k] = dq.get(k, 0) + 1
    print("%s: qlen(orig) - qlen(nostar) over %d queries with hits in both runs: %s | lengths without '*': %s"
          % (lab, len(both), dict(sorted(dq.items())), "qlen(nostar) == length for all" if all(ql["nostar"][x] == L[x] for x in ql["nostar"]) else "qlen(nostar) != length for some"))
bo, bn_ = B[("EU15580", "orig")]["either"], B[("EU15580", "nostar")]["either"]
Ho = {x for x, e in bo.items() if e <= EMAX}; Hn = {x for x, e in bn_.items() if e <= EMAX}
if len(Ho) != 218: raise SystemExit("ABORT: original EU15580 HIT is %d, expected 218" % len(Ho))
L15 = lengths(R + "prostt5_main_eu15580_nostar.faa")
with open(O3, "w") as o:
    o.write("rep\tlength\tbestE_orig\tbestE_nostar\tstatus\n")
    for x in sorted(Ho | Hn, key=lambda x: (x not in Ho, x not in Hn, x)):
        st = "kept" if x in Ho and x in Hn else ("lost" if x in Ho else "gained")
        o.write("%s\t%d\t%s\t%s\t%s\n" % (x, L15[x], "%.3g" % bo[x] if x in bo else "none", "%.3g" % bn_[x] if x in bn_ else "none", st))
for p, t in ((O1, t1), (O2, t2)):
    with open(p, "w") as fh:
        for r in t: fh.write("\t".join(map(str, r)) + "\n")
for t in (t1, t2):
    for r in t: print("\t".join(map(str, r)))
    print()
print("EU15580 HIT (either, E <= 1e-3): original %d, nostar %d | kept %d, lost %d, gained %d"
      % (len(Ho), len(Hn), len(Ho & Hn), len(Ho - Hn), len(Hn - Ho)))
print("wrote", O1); print("wrote", O2); print("wrote", O3)
