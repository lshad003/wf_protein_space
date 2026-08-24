# Step 16b: materialize the support-filtered cluster set at the 95% tier.
# Criterion: >= 3 members AND members from >= 2 metagenomes.
# Cross-checks: 40,550,595 member lines, 6,182,117 clusters, 2,069,453 supported.
import sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
IN  = W + "/results/clusters_95.tsv"
OUT = W + "/results/supported_reps_95.txt"
SUM = W + "/results/supported_reps_95.summary.txt"
MULTI = "\x00MULTI"
cnt = {}
first = {}
n = 0
with open(IN) as fh:
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        rep, mem = p[0], p[1]
        stem = mem.split("__", 1)[0]
        cnt[rep] = cnt.get(rep, 0) + 1
        f = first.get(rep)
        if f is None:
            first[rep] = sys.intern(stem)
        elif f != MULTI and f != stem:
            first[rep] = MULTI
        n += 1
        if n % 5000000 == 0:
            print("lines:", n, flush=True)
supported = sorted(r for r in cnt if cnt[r] >= 3 and first[r] == MULTI)
with open(OUT, "w") as out:
    out.write("\n".join(supported) + "\n")
lines = ["criterion: >=3 members AND >=2 metagenomes, 95%% tier",
         "member_lines\t%d\texpect\t40550595" % n,
         "clusters\t%d\texpect\t6182117" % len(cnt),
         "supported\t%d\texpect\t2069453" % len(supported)]
with open(SUM, "w") as out:
    out.write("\n".join(lines) + "\n")
print("\n".join(lines))
ok = (n == 40550595 and len(cnt) == 6182117 and len(supported) == 2069453)
print("CROSS-CHECK:", "PASS" if ok else "MISMATCH, stop and investigate")
