# Step 34c: strict EU 30% families, recurrence, NMPFams2 and the Cse2 positive control.
# Usage: python3 step34c_families.py <clusters.tsv from mmseqs createtsv (family rep, member rep)>
# Inputs: results/step34a_rep_origin.tsv (strict EU = EU and Tiara-prokaryotic), results/clusters_95_n88.tsv
#   (95% rep, member protein), metadata/samples_88_treatment.tsv (stem, cohort), results/unk_filtered_aplakidou.faa
#   (step26 NMPFams2 query set), results/nmpfams2_{metag,iso}_rep_hits.tsv (matched = any row, as step26d4/step32b).
# Per family: member reps (= distinct 95% variants; each catalog rep is one 95% cluster), member proteins of those
#   reps from the 88 samples, metagenomes, animals (stem prefix before ".", as step26f2), cohorts.
# Recurrent = >= 10 animals AND >= 2 cohorts (step26f2 criterion). Single-variant = family with one member rep.
# Outputs: results/step34c_strict_eu_families.tsv, results/step34c_strict_eu_nmpfams2.tsv; refuses to overwrite.
import os, sys, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results/"
OUT1 = R + "step34c_strict_eu_families.tsv"
OUT2 = R + "step34c_strict_eu_nmpfams2.tsv"
for f in (OUT1, OUT2):
    if os.path.exists(f): raise SystemExit("REFUSING: " + f + " exists")
CSE2 = "UHM636.41030__k141_332084_5"
eu = set()
with open(R + "step34a_rep_origin.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[1] == "EU" and f[2] == "Tiara-prokaryotic": eu.add(f[0])
print("strict EU reps:", len(eu), flush=True)
if len(eu) != 41392: raise SystemExit("ABORT: expected 41,392")
fam_of = {}; fam = collections.defaultdict(list)
for l in open(sys.argv[1]):
    a, b = l.rstrip("\n").split("\t")
    if b in fam_of: raise SystemExit("ABORT: rep in two families: " + b)
    fam_of[b] = a; fam[a].append(b)
print("clustered reps:", len(fam_of), "| families:", len(fam), flush=True)
if set(fam_of) != eu: raise SystemExit("ABORT: clustered reps differ from strict EU set")
coh = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); coh[f[0]] = f[1]
print("stems:", len(coh), dict(collections.Counter(coh.values())), flush=True)
mem = collections.defaultdict(list); off = 0; rep_of = {}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m == CSE2: rep_of[CSE2] = r
        if r not in eu: continue
        s = m.split("__", 1)[0]
        if s not in coh: off += 1; continue
        mem[r].append(s)
print("member proteins of strict EU reps from the 88 samples:", sum(len(v) for v in mem.values()),
      "| members with a stem outside the 88:", off, "| reps without a member:", len(eu - set(mem)), flush=True)
rows = []
for f, rs in fam.items():
    st = [s for r in rs for s in mem.get(r, [])]
    S = set(st)
    an = {s.split(".", 1)[0] for s in S}; co = {coh[s] for s in S}
    rec = len(an) >= 10 and len(co) >= 2
    rows.append((f, len(rs), len(rs), len(st), len(S), len(an), len(co), ",".join(sorted(co)), int(rec)))
rows.sort(key=lambda x: (-x[8], -x[5], -x[3], x[0]))
with open(OUT1, "w") as fh:
    fh.write("family_rep\tn_member_reps\tn_variants95\tn_member_proteins\tn_metagenomes\tn_animals\tn_cohorts\tcohorts\trecurrent\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
N = len(rows); rec = [r for r in rows if r[8]]
sv = sum(1 for r in rec if r[2] == 1); sva = sum(1 for r in rows if r[2] == 1)
print("\nfamilies: %d | with >= 2 member reps: %d | single-variant (all families): %d (%.2f%%)" % (N, N - sva, sva, 100.0 * sva / N))
print("recurrent (>= 10 animals AND >= 2 cohorts): %d of %d (%.2f%%)" % (len(rec), N, 100.0 * len(rec) / N))
if rec:
    print("single-variant among recurrent: %d of %d (%.2f%%) [old 15,580 set: 14,918, 95.75%%]" % (sv, len(rec), 100.0 * sv / len(rec)))
    print("strict EU reps in recurrent families: %d of %d" % (sum(r[1] for r in rec), len(eu)))
for lab, fn in ((">= 10 animals", lambda r: r[5] >= 10), (">= 2 cohorts", lambda r: r[6] >= 2),
                ("all 3 cohorts", lambda r: r[6] == 3), (">= 10 metagenomes AND >= 2 cohorts", lambda r: r[4] >= 10 and r[6] >= 2)):
    print("  %s: %d" % (lab, sum(1 for r in rows if fn(r))))
print("member reps per family (all): " + ", ".join("%s %d" % (k, sum(1 for r in rows if lo <= r[1] < hi))
      for k, lo, hi in (("1", 1, 2), ("2", 2, 3), ("3-4", 3, 5), ("5-9", 5, 10), (">=10", 10, 10**9))))
if rec:
    print("member reps per recurrent family: " + ", ".join("%s %d" % (k, sum(1 for r in rec if lo <= r[1] < hi))
          for k, lo, hi in (("1", 1, 2), ("2", 2, 3), ("3-4", 3, 5), ("5-9", 5, 10), (">=10", 10, 10**9))))
# NMPFams2
q = {l[1:].split()[0] for l in open(R + "unk_filtered_aplakidou.faa") if l[0] == ">"}
hm = {l.split("\t", 1)[0] for l in open(R + "nmpfams2_metag_rep_hits.tsv")}
hi = {l.split("\t", 1)[0] for l in open(R + "nmpfams2_iso_rep_hits.tsv")}
recf = {r[0] for r in rec}
with open(OUT2, "w") as fh:
    fh.write("rep\tfamily_rep\trecurrent_family\tin_step26_set\tnmp_metag\tnmp_iso\tmatched\n")
    for r in sorted(eu):
        i = r in q
        fh.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\n" % (r, fam_of[r], int(fam_of[r] in recf), int(i),
                 int(r in hm) if i else "NA", int(r in hi) if i else "NA", int(r in hm or r in hi) if i else "NA"))
print("\nNMPFams2 (step26 query set results/unk_filtered_aplakidou.faa, %d; matched = any metag or iso rep hit)" % len(q))
for lab, S in (("all strict EU reps", eu), ("reps in recurrent families", {r for r in eu if fam_of[r] in recf}),
               ("reps in non-recurrent families", {r for r in eu if fam_of[r] not in recf})):
    I = S & q; u = sum(1 for r in I if r not in hm and r not in hi)
    print("  %-32s n=%6d | in step26 set %6d (%.2f%%) | unmatched %6d (%s) | metag hit %d, iso hit %d"
          % (lab, len(S), len(I), 100.0 * len(I) / len(S) if S else 0, u, "%.2f%%" % (100.0 * u / len(I)) if I else "NA",
             len(I & hm), len(I & hi)))
print("  [old: prokaryotic >= 3 kb EU 23,781/24,415 unmatched (97.40%); all EU 265,061/266,848 (99.33%)]")
fu = {}
for r in eu:
    if r in q: fu.setdefault(fam_of[r], []).append(r not in hm and r not in hi)
nf = len(fu); fm = sum(1 for v in fu.values() if not all(v))
print("  families with >= 1 queried rep: %d | with any matched rep: %d | all queried reps unmatched: %d (%.2f%%)" % (nf, fm, nf - fm, 100.0 * (nf - fm) / nf))
# positive control
print("\npositive control", CSE2)
r = rep_of.get(CSE2)
print("  95% rep:", r, "| is itself a rep:", CSE2 == r, "| rep strict EU:", r in eu)
if r in eu:
    f = fam_of[r]; x = [y for y in rows if y[0] == f][0]
    print("  family %s: member reps %d, member proteins %d, metagenomes %d, animals %d, cohorts %s, recurrent %d"
          % (f, x[1], x[3], x[4], x[5], x[7], x[8]))
    print("  rep in step26 set:", r in q, "| NMPFams2 matched:", (r in hm or r in hi) if r in q else "NA")
    print("  PASS" if x[8] else "  FAIL: not recurrent")
else:
    print("  FAIL: rep not in strict EU set")
print("wrote", OUT1, "rows:", N)
print("wrote", OUT2, "rows:", len(eu))
