# Step 36b prep: search set for the Pfam near-miss scan of the 4,509 recurrent EU families.
# Members: every member of the families' 95% clusters (results/step34h/strict_members.tsv, strict 0 or 1; family from
#   results/step34c_strict_eu_nmpfams2.tsv; recurrent == 1 and excluded == 0 in step34c2c_strict_v2_families.tsv),
#   88 metagenomes only (metadata/samples_88_treatment.tsv). Strict members are flagged; non-strict members are searched
#   too so the all-member counts (e.g. Cse2, 45 members) can be reported beside the strict ones. Family reps that are
#   not in the 88 (UHM586.41010) are searched as the representative only (strict = rep_not_in_88).
# Sequences: samtools faidx on catalog/db/LsPS_AA_all.fasta.gz (run by step36b_prep.sh), terminal "*" removed;
#   identical sequences searched once (results/step36b/uniq_map.tsv: uid -> member).
# Outputs (refuse to overwrite): results/step36b/members.tsv, ids.txt, then (stage split) uniq.faa, uniq_map.tsv,
#   chunks/C.<k>.faa.
import os, sys, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step36b/"
STAGE = sys.argv[1]
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
if STAGE == "members":
    for p in (D + "members.tsv", D + "ids.txt"):
        if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
    S88 = {r["stem"] for r in tsv(W + "/metadata/samples_88_treatment.tsv")}
    fams = {r["family_rep"] for r in tsv(R + "step34c2c_strict_v2_families.tsv") if r["recurrent"] == "1" and r["excluded"] == "0"}
    r2f = {r["rep"]: r["family_rep"] for r in tsv(R + "step34c_strict_eu_nmpfams2.tsv") if r["family_rep"] in fams}
    if len(fams) != 4509 or len(S88) != 88: raise SystemExit("ABORT: family or sample count")
    M = []; drop = 0
    for r in tsv(R + "step34h/strict_members.tsv"):
        if r["rep"] not in r2f: continue
        if r["member"].split("__", 1)[0] not in S88: drop += 1; continue
        M.append((r2f[r["rep"]], r["rep"], r["member"], r["strict"]))
    have = {m[2] for m in M}
    nofr = [f for f in fams if f not in have]
    # 2026-10-23: 26 family reps are UHM586.41010 proteins (reps chosen with 89 samples); the rep sequence is still
    #   searched (strict = "rep_not_in_88") but never counted as a member.
    if any(not f.startswith("UHM586.") for f in nofr): raise SystemExit("ABORT: family reps missing for another reason")
    for f in nofr: M.append((f, f, f, "rep_not_in_88")); have.add(f)
    print("family reps outside the 88 (searched, not counted as members): %d" % len(nofr))
    with open(D + "members.tsv", "w") as o:
        o.write("family_rep\trep95\tmember\tstrict\n")
        for x in sorted(M): o.write("\t".join(x) + "\n")
    with open(D + "ids.txt", "w") as o:
        for m in sorted(have): o.write(m + "\n")
    c = collections.Counter(x[3] for x in M)
    print("families %d | 95%% reps %d | rows %d (strict %d, non-strict %d) | dropped (not in the 88) %d | distinct ids %d"
          % (len(fams), len(r2f), len(M), c["1"], c["0"], drop, len(have)))
elif STAGE == "split":
    NCH = int(sys.argv[2])
    for p in (D + "uniq.faa", D + "uniq_map.tsv", D + "chunks"):
        if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
    seq = {}; h = None; buf = []
    for l in open(D + "members.faa"):
        if l[0] == ">":
            if h: seq[h] = "".join(buf).rstrip("*")
            h = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if h: seq[h] = "".join(buf).rstrip("*")
    ids = [l.strip() for l in open(D + "ids.txt")]
    miss = [i for i in ids if i not in seq]
    if miss: raise SystemExit("ABORT: %d ids not extracted, e.g. %s" % (len(miss), miss[:3]))
    u = {}
    for i in ids: u.setdefault(seq[i], []).append(i)
    U = sorted(u.items(), key=lambda t: (-len(t[0]), t[1][0]))
    os.makedirs(D + "chunks")
    outs = [open(D + "chunks/C.%d.faa" % k, "w") for k in range(NCH)]
    with open(D + "uniq.faa", "w") as of, open(D + "uniq_map.tsv", "w") as om:
        om.write("uid\tmember\n")
        for k, (s, ms) in enumerate(U):
            uid = "S%06d" % k
            of.write(">%s\n%s\n" % (uid, s)); outs[k % NCH].write(">%s\n%s\n" % (uid, s))
            for m in ms: om.write("%s\t%s\n" % (uid, m))
    for o in outs: o.close()
    print("members %d | distinct sequences %d | residues %d | chunks %d" % (len(ids), len(U), sum(len(s) for s, _ in U), NCH))
