# Step 34c2c part 2: strict_v2 = strict (Step 34a: Tiara-prokaryotic, not fly95 or frog95) minus reps whose contig is
#   (a) fly at >= 90% identity with the same coverage rule (minimap2 -x asm10 --secondary=no -c to dm6; merged query
#       coverage >= 0.50 AND sum(matches)/sum(block length) >= 0.90; PAFs results/step33g_paf for step33c_lists contigs,
#       results/step34a_paf for step34a_lists contigs), or
#   (b) Eukaryota by the v2 contig LCA domain (results/gene_taxonomy_n88_v2.tsv; 48 of 88 stems have contig taxonomy).
# Then the Step 34c family statistics on strict_v2 by dropping removed reps from the existing 30% families
#   (results/step34c_strict_eu_clusters30.tsv; no re-clustering): members, animals, cohorts recomputed from kept reps.
# Outputs: results/step34a_rep_origin_v2.tsv (step34a columns + fly90, lca_domain_v2, strict_v2, removed_by; old file
#   kept); rows appended to results/step34c2c_summary.tsv (written by step34c2c_taxv2.py).
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = R + "step34a_rep_origin_v2.tsv"; SUM = R + "step34c2c_summary.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
if not os.path.exists(SUM) or "step33e" not in open(SUM).read(): raise SystemExit("ABORT: run step34c2c_taxv2.py first")
if "strict_v2" in open(SUM).read(): raise SystemExit("REFUSING: part 2 rows already in " + SUM)
summ = []
def S(sec, metric, old, new): summ.append((sec, metric, str(old), str(new))); print("  %-62s old %-24s new %s" % (metric, old, new), flush=True)
def contig(g):
    s, rest = g.split("__", 1); return (s, rest.rsplit("_", 1)[0])
rows = []; hdr = None
with open(R + "step34a_rep_origin.tsv") as fh:
    hdr = next(fh).rstrip("\n").split("\t")
    for l in fh: rows.append(l.rstrip("\n").split("\t"))
print("step34a rows %d | columns %s" % (len(rows), hdr), flush=True)
# fly calls at 95 and 90 from the PAFs
lists = collections.defaultdict(set)
stems = [l.strip() for l in open(R + "step34a_lists/stems.txt") if l.strip()]
def calls(p, s, want):
    iv, mt, bl, ql = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
    if os.path.exists(p):
        for l in open(p):
            f = l.split("\t"); k = (s, f[0])
            if k not in want: continue
            ql[k] = int(f[1]); iv[k].append((int(f[2]), int(f[3]))); mt[k] += int(f[9]); bl[k] += int(f[10])
    out = {}
    for k in want:
        if k not in iv: out[k] = (0.0, 0.0); continue
        v = sorted(iv[k]); c, cs, ce = 0, -1, -1
        for a, b in v:
            if a > ce: c += max(0, ce - cs); cs, ce = a, b
            else: ce = max(ce, b)
        out[k] = ((c + ce - cs) / ql[k], mt[k] / bl[k])
    return out
cv = {}
for s in stems:
    a = {(s, l.split("\t")[0]) for l in open(R + "step33c_lists/%s.tsv" % s) if l.strip()} if os.path.exists(R + "step33c_lists/%s.tsv" % s) else set()
    b = {(s, l.strip()) for l in open(R + "step34a_lists/%s.txt" % s) if l.strip()}
    if a & b: raise SystemExit("ABORT: contig in both lists " + s)
    cv.update(calls(R + "step33g_paf/%s.fly.paf" % s, s, a)); cv.update(calls(R + "step34a_paf/%s.fly.paf" % s, s, b))
f95 = lambda k: cv[k][0] >= 0.5 and cv[k][1] >= 0.95
f90 = lambda k: cv[k][0] >= 0.5 and cv[k][1] >= 0.90
miss = sum(1 for r in rows if contig(r[0]) not in cv)
bad95 = sum(1 for r in rows if contig(r[0]) in cv and int(f95(contig(r[0]))) != int(r[4]))
g33 = {}
for i, l in enumerate(open(R + "step33g_contig_fly.tsv")):
    if i: f = l.rstrip("\n").split("\t"); g33[(f[0], f[1])] = f[7] == "1"
bad90 = sum(1 for k, v in g33.items() if k in cv and f90(k) != v)
print("rep contigs without a fly call %d | fly95 disagreements vs step34a %d | fly90 disagreements vs step33g_contig_fly %d" % (miss, bad95, bad90), flush=True)
if miss or bad95 or bad90: raise SystemExit("ABORT: fly calls do not reproduce")
dv2 = {}
with open(R + "gene_taxonomy_n88_v2.tsv") as fh:
    next(fh)
    for l in fh: p = l.rstrip("\n").split("\t"); dv2[p[0]] = p[4]
# strict_v2
rm = collections.defaultdict(collections.Counter); st = collections.Counter(); keep = set(); lab = {}
with open(OUT, "w") as o:
    o.write("\t".join(hdr + ["fly90", "lca_domain_v2", "strict_v2", "removed_by"]) + "\n")
    for r in rows:
        k = contig(r[0]); a = f90(k); d = dv2.get(r[0], "NA"); strict = r[2] == "Tiara-prokaryotic"
        why = "-"
        if strict:
            st[r[1]] += 1; e = d == "Eukaryota"
            why = "fly90+euk_lca" if a and e else ("fly90" if a else ("euk_lca" if e else "kept"))
            rm[r[1]][why] += 1
            if why == "kept": keep.add(r[0])
            lab[r[0]] = why
        o.write("\t".join(r + [str(int(a)), d, str(int(strict and why == "kept")), why]) + "\n")
print("\nstrict -> strict_v2 by class (removed: fly90 only / euk_lca only / both)")
for c in ("K", "KWP", "GU", "EU", "all"):
    n = sum(st.values()) if c == "all" else st[c]
    x = lambda w: sum(rm[k][w] for k in rm) if c == "all" else rm[c][w]
    S("strict_v2", "%s strict %d: removed fly90 only / euk_lca only / both" % (c, n), "-",
      "%d / %d / %d" % (x("fly90"), x("euk_lca"), x("fly90+euk_lca")))
    S("strict_v2", "%s reps" % c, n, "%d (%.2f%% removed)" % (x("kept"), 100.0 * (n - x("kept")) / n))
tot = sum(st.values()); kv = len(keep)
S("strict_v2", "EU share of set", "%.2f%%" % (100.0 * st["EU"] / tot), "%.2f%%" % (100.0 * rm["EU"]["kept"] / kv))
# families
fam_of = {}; reps_of = collections.defaultdict(list)
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t"); fam_of[b] = a; reps_of[a].append(b)
EU1 = set(fam_of)
if len(EU1) != 41392: raise SystemExit("ABORT: step34c families")
EU2 = EU1 & keep
print("strict EU reps %d -> strict_v2 EU %d" % (len(EU1), len(EU2)), flush=True)
coh = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: f = l.rstrip("\n").split("\t"); coh[f[0]] = f[1]
mem = collections.defaultdict(list); rep_of = {}
CSE2 = "UHM636.41030__k141_332084_5"
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m == CSE2: rep_of[m] = r
        if r in EU1: mem[r].append(m.split("__", 1)[0])
def stats(K):
    out = {}
    for f, rs in reps_of.items():
        rs = [r for r in rs if r in K]
        if not rs: continue
        S_ = {s for r in rs for s in mem[r]}
        an = {s.split(".", 1)[0] for s in S_}; co = {coh[s] for s in S_}
        out[f] = (len(rs), len(an), len(co), len(an) >= 10 and len(co) >= 2)
    return out
A, B = stats(EU1), stats(EU2)
ra = {f for f, v in A.items() if v[3]}; rb = {f for f, v in B.items() if v[3]}
S("34c families", "families", len(A), len(B))
S("34c families", "single-rep families (all)", "%d (%.2f%%)" % (sum(v[0] == 1 for v in A.values()), 100.0 * sum(v[0] == 1 for v in A.values()) / len(A)),
  "%d (%.2f%%)" % (sum(v[0] == 1 for v in B.values()), 100.0 * sum(v[0] == 1 for v in B.values()) / len(B)))
S("34c families", "recurrent (>= 10 animals AND >= 2 cohorts)", len(ra), len(rb))
S("34c families", "reps in recurrent families", sum(A[f][0] for f in ra), sum(B[f][0] for f in rb))
S("34c families", "single-variant among recurrent", "%d (%.2f%%)" % (sum(A[f][0] == 1 for f in ra), 100.0 * sum(A[f][0] == 1 for f in ra) / len(ra)),
  "%d (%.2f%%)" % (sum(B[f][0] == 1 for f in rb), 100.0 * sum(B[f][0] == 1 for f in rb) / len(rb)))
S("34c families", "recurrent families lost entirely / lost recurrence only", "-",
  "%d / %d" % (sum(1 for f in ra if f not in B), sum(1 for f in ra if f in B and f not in rb)))
q = {l[1:].split()[0] for l in open(R + "unk_filtered_aplakidou.faa") if l[0] == ">"}
hit = {l.split("\t", 1)[0] for x in ("metag", "iso") for l in open(R + "nmpfams2_%s_rep_hits.tsv" % x)}
def nmp(E, rec, fam):
    I = E & q if fam is None else {r for r in E & q if fam_of[r] in rec}
    u = sum(1 for r in I if r not in hit)
    return "%d/%d (%.2f%%)" % (u, len(I), 100.0 * u / len(I))
S("34c NMPFams2", "unmatched, strict EU reps searched", nmp(EU1, None, None), nmp(EU2, None, None))
S("34c NMPFams2", "unmatched, reps in recurrent families", nmp(EU1, ra, 1), nmp(EU2, rb, 1))
dro = []
with open(R + "step34c2_recurrent_genus.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[0] == "recurrent_EU_family" and f[7] == "Drosophila": dro.append(f[0 + 1])
gone = sum(1 for f in dro if f not in B); norec = sum(1 for f in dro if f in B and f not in rb)
why = collections.Counter(lab[r] for f in dro for r in reps_of[f])
S("34c Drosophila-majority", "22 recurrent families: gone / present not recurrent / still recurrent", "0 / 0 / %d" % len(dro),
  "%d / %d / %d" % (gone, norec, len(dro) - gone - norec))
S("34c Drosophila-majority", "their reps by rule", "-", ";".join("%s:%d" % kv for kv in sorted(why.items())))
for f in dro:
    if f in rb: print("   still recurrent:", f, [(r, lab[r]) for r in reps_of[f]])
r = rep_of[CSE2]; f = fam_of.get(r)
ok = r in EU2 and f in rb
S("34c control", "Cse2 rep %s: strict_v2, family recurrent (members/animals/cohorts)" % r,
  "kept, recurrent (%d reps, %d animals, %d cohorts)" % A[f][:3], "%s, %s%s" % ("kept" if r in EU2 else "REMOVED " + lab[r],
  "recurrent" if f in rb else "NOT recurrent", " (%d reps, %d animals, %d cohorts)" % B[f][:3] if f in B else ""))
print("  Cse2:", "PASS" if ok else "FAIL")
with open(SUM, "a") as o:
    for x in summ: o.write("\t".join(x) + "\n")
print("wrote", OUT, "| appended part 2 to", SUM)
