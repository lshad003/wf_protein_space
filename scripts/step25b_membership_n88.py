# Step 25b: cluster membership at 88 metagenomes. Drops every member row from
# UHM586.41010 at each tier and drops clusters left empty. No reclustering;
# representative IDs are kept as cluster labels even where the representative
# sequence itself came from UHM586.41010. Also writes two n88 intermediates:
#   results/reps_95_n88.txt          95% clusters that keep at least one member
#   results/unk_clusters_50_n88.tsv  unk_clusters_50.tsv restricted to members
#                                    in results/supported_reps_95_n88.txt
import collections, sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space/results/"
PFX = "UHM586.41010__"
res = {}
for tier in ("95", "90", "50"):
    n_in = collections.Counter(); n_out = collections.Counter()
    rin = rout = 0
    with open(W + "clusters_%s.tsv" % tier) as fi, open(W + "clusters_%s_n88.tsv" % tier, "w") as fo:
        for line in fi:
            rep, mem = line.split("\t", 2)[:2]
            rin += 1; n_in[rep] += 1
            if mem.startswith(PFX): continue
            fo.write(line); rout += 1; n_out[rep] += 1
    cin, cout = len(n_in), len(n_out)
    rep_from_586 = sum(1 for r in n_out if r.startswith(PFX))
    print("tier %s: rows in %d out %d dropped %d | clusters in %d out %d dropped %d | kept clusters whose rep is a UHM586.41010 protein %d"
          % (tier, rin, rout, rin - rout, cin, cout, cin - cout, rep_from_586), flush=True)
    res[tier] = (rin - rout, cin - cout)
    if tier == "95":
        with open(W + "reps_95_n88.txt", "w") as o:
            for r in sorted(n_out): o.write(r + "\n")
        print("wrote reps_95_n88.txt:", cout)
    del n_in, n_out
if res["95"] != (292633, 10515):
    sys.exit("ABORT: 95%% tier lost %d rows and %d clusters, expected 292633 and 10515" % res["95"])
print("CHECK 95% tier: PASS")
sup = set(open(W + "supported_reps_95_n88.txt").read().split())
fam_in = collections.Counter(); fam_out = collections.Counter(); mi = mo = 0
with open(W + "unk_clusters_50.tsv") as fi, open(W + "unk_clusters_50_n88.tsv", "w") as fo:
    for line in fi:
        r, m = line.rstrip("\n").split("\t")[:2]
        mi += 1; fam_in[r] += 1
        if m in sup:
            fo.write(line); mo += 1; fam_out[r] += 1
ge3 = lambda c: sum(1 for v in c.values() if v >= 3)
print("unk_clusters_50: members in %d out %d | families in %d out %d | families >=3 members in %d out %d"
      % (mi, mo, len(fam_in), len(fam_out), ge3(fam_in), ge3(fam_out)))
print("DONE step25b")
