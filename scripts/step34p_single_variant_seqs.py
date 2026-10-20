# Step 34p: distinct member sequences in the 4,041 single-variant recurrent families (recurrent 1, excluded 0,
#   n_member_reps 1 in results/step34c2c_strict_v2_families.tsv; "single-variant" = one 95% representative,
#   scripts/step34c2c_exclude.py:67-68). Family -> rep from results/step34c_strict_eu_clusters30.tsv; members, their
#   distinct-sequence IDs (uid) and strict eligibility from results/step34h/strict_members.tsv (rep, member, uid, strict).
# Output: results/step34p_single_variant_seqs.tsv (per family) ; refuses to overwrite.
import os, collections
R = "/bigdata/stajichlab/lshad003/wf_protein_space/results/"
OUT = R + "step34p_single_variant_seqs.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
single = set()
for l in list(open(R + "step34c2c_strict_v2_families.tsv"))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[8] == "1" and f[9] == "0" and f[1] == "1": single.add(f[0])
if len(single) != 4041: raise SystemExit("ABORT: %d single-variant families" % len(single))
rep_of = {}
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t")
    if a in single: rep_of[a] = b
if len(rep_of) != 4041: raise SystemExit("ABORT: %d family reps mapped" % len(rep_of))
reps = set(rep_of.values())
nm = collections.Counter(); ns = collections.Counter(); ua = collections.defaultdict(set); us = collections.defaultdict(set)
for i, l in enumerate(open(R + "step34h/strict_members.tsv")):
    if i == 0: continue
    r, m, u, s = l.rstrip("\n").split("\t")
    if r not in reps: continue
    nm[r] += 1; ua[r].add(u)
    if s == "1": ns[r] += 1; us[r].add(u)
if len(ua) != 4041: raise SystemExit("ABORT: %d reps with members" % len(ua))
with open(OUT, "w") as o:
    o.write("family_rep\trep\tn_members\tn_strict_members\tdistinct_seqs_all\tdistinct_seqs_strict\n")
    for f in sorted(rep_of):
        r = rep_of[f]; o.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (f, r, nm[r], ns[r], len(ua[r]), len(us[r])))
A = [len(ua[r]) for r in reps]; S = [len(us[r]) for r in reps]
print("single-variant recurrent families: %d" % len(reps))
print("distinct sequences, all members: > 1 in %d, exactly 1 in %d, median %d, max %d" % (sum(x > 1 for x in A), sum(x == 1 for x in A), sorted(A)[len(A) // 2], max(A)))
print("distinct sequences, strict-eligible members: > 1 in %d, exactly 1 in %d, 0 in %d" % (sum(x > 1 for x in S), sum(x == 1 for x in S), sum(x == 0 for x in S)))
print("wrote", OUT)
