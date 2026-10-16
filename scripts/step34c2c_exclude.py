# Step 34c2c exclusion: the strict_v2 family table, with the 5 Drosophila-majority recurrent families that survived
# strict_v2 (step34c2b hard cases) flagged as excluded.
# strict_v2 reps: results/step34a_rep_origin_v2.tsv (strict_v2 = 1, class4 EU). Families: results/step34c_strict_eu_clusters30.tsv
#   with removed reps dropped (no re-clustering), statistics as step34c_families.py (members via clusters_95_n88.tsv,
#   animal = stem prefix before ".", cohorts from metadata/samples_88_treatment.tsv; recurrent = >= 10 animals AND >= 2
#   cohorts; NMPFams2 matched = any metag or iso rep hit). results/step34a_rep_origin_v2.tsv is not changed.
# Rep level: results/step34a_rep_origin_v2_final.tsv = step34a_rep_origin_v2.tsv + strict_final (strict_v2 minus the reps
#   of the excluded families). All later steps use strict_final.
# Outputs (refuse to overwrite): results/step34c2c_excluded_families.tsv, results/step34c2c_strict_v2_families.tsv,
#   results/step34a_rep_origin_v2_final.tsv
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT1 = R + "step34c2c_excluded_families.tsv"; OUT2 = R + "step34c2c_strict_v2_families.tsv"
OUT3 = R + "step34a_rep_origin_v2_final.tsv"
for p in (OUT1, OUT2, OUT3):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
EXCL = ["UHM38.35752__k141_20484_2", "UHM38.35752__k141_20484_6", "UHM558.40992__k141_126719_8",
        "UHM626.41025__k141_36605_3", "UHM633.41028__k141_435107_1"]
WHY = "Drosophila majority genus; failed fly mapping rule"
CSE2 = "UHM636.41030__k141_332084_5"
keep = set()
with open(R + "step34a_rep_origin_v2.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, isv = h.index("class4"), h.index("strict_v2")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[ic] == "EU" and f[isv] == "1": keep.add(f[0])
print("strict_v2 EU reps:", len(keep), flush=True)
if len(keep) != 41119: raise SystemExit("ABORT: expected 41,119")
reps_of = collections.defaultdict(list); fam_of = {}
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t")
    if b in keep: reps_of[a].append(b); fam_of[b] = a
if set(fam_of) != keep: raise SystemExit("ABORT: strict_v2 EU reps missing from step34c clusters")
dro = set()
with open(R + "step34c2_recurrent_genus.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[0] == "recurrent_EU_family" and f[7] == "Drosophila": dro.add(f[1])
if not set(EXCL) <= dro or not set(EXCL) <= set(reps_of): raise SystemExit("ABORT: excluded families are not surviving Drosophila-majority families")
if sorted(f for f in dro if f in reps_of) != sorted(EXCL): raise SystemExit("ABORT: surviving Drosophila-majority families differ from the list")
coh = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: f = l.rstrip("\n").split("\t"); coh[f[0]] = f[1]
mem = collections.defaultdict(list); rep_of = {}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m == CSE2: rep_of[m] = r
        if r in keep: mem[r].append(m.split("__", 1)[0])
rows = []
for f, rs in reps_of.items():
    st = [s for r in rs for s in mem[r]]; S = set(st)
    an = {s.split(".", 1)[0] for s in S}; co = {coh[s] for s in S}
    rec = len(an) >= 10 and len(co) >= 2; ex = f in EXCL
    rows.append((f, len(rs), len(rs), len(st), len(S), len(an), len(co), ",".join(sorted(co)), int(rec), int(ex), WHY if ex else "-"))
rows.sort(key=lambda x: (-x[8], -x[5], -x[3], x[0]))
with open(OUT2, "w") as o:
    o.write("family_rep\tn_member_reps\tn_variants95\tn_member_proteins\tn_metagenomes\tn_animals\tn_cohorts\tcohorts\trecurrent\texcluded\texclusion_reason\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
with open(OUT1, "w") as o:
    o.write("family_rep\treps\treason\n")
    for f in EXCL: o.write("%s\t%s\t%s\n" % (f, ",".join(sorted(reps_of[f])), WHY))
ok = [r for r in rows if not r[9]]; rec = [r for r in ok if r[8]]
print("families %d | excluded %d | final families %d" % (len(rows), len(EXCL), len(ok)))
print("recurrent before exclusion %d | final recurrent %d (reps %d) | single-variant %d (%.2f%%)" % (
    sum(r[8] for r in rows), len(rec), sum(r[1] for r in rec), sum(r[1] == 1 for r in rec), 100.0 * sum(r[1] == 1 for r in rec) / len(rec)))
q = {l[1:].split()[0] for l in open(R + "unk_filtered_aplakidou.faa") if l[0] == ">"}
hit = {l.split("\t", 1)[0] for x in ("metag", "iso") for l in open(R + "nmpfams2_%s_rep_hits.tsv" % x)}
for lab, F in (("all non-excluded", {r[0] for r in ok}), ("recurrent non-excluded", {r[0] for r in rec})):
    I = {x for f in F for x in reps_of[f]} & q; u = sum(1 for x in I if x not in hit)
    print("NMPFams2 unmatched, %s: %d/%d (%.2f%%)" % (lab, u, len(I), 100.0 * u / len(I)))
r = rep_of[CSE2]; f = fam_of.get(r); x = [y for y in rows if y[0] == f]
print("Cse2 rep %s family %s:" % (r, f), x[0] if x else None)
print("Cse2:", "PASS" if x and x[0][8] and not x[0][9] else "FAIL")
xr = {x for f in EXCL for x in reps_of[f]}
n = collections.Counter(); nv2 = collections.Counter()
with open(R + "step34a_rep_origin_v2.tsv") as fh, open(OUT3, "w") as o:
    o.write(next(fh).rstrip("\n") + "\tstrict_final\n")
    for l in fh:
        f = l.rstrip("\n").split("\t"); v2 = f[isv] == "1"; fin = v2 and f[0] not in xr
        nv2[f[ic]] += v2; n[f[ic]] += fin
        o.write(l.rstrip("\n") + "\t%d\n" % fin)
print("\nrep level: excluded reps %d" % len(xr))
for c in ("K", "KWP", "GU", "EU"):
    print("  %-4s strict_v2 %9d -> strict_final %9d (%.2f%%)" % (c, nv2[c], n[c], 100.0 * n[c] / sum(n.values())))
print("  all  strict_v2 %9d -> strict_final %9d" % (sum(nv2.values()), sum(n.values())))
if n["EU"] != 41114: raise SystemExit("ABORT: strict_final EU %d, expected 41,114" % n["EU"])
print("wrote", OUT1, OUT2, "and", OUT3)
