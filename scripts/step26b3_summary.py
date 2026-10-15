# Step 26b3 summary: represented vs not represented in NMPFams2.
# Run with /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 after
# scripts/step26b3_nmpfams2_search.sh; covers every source whose hit files exist.
# Hit criterion: mmseqs search E <= 1e-3, >=30% identity, >=70% coverage of both
# query and target (--cov-mode 0). A family counts as represented when any of its
# queried 95% representatives has a hit.
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results"
OUT = R + "/nmpfams2_summary.tsv"

queries = [l[1:].split()[0] for l in open(R + "/unk_filtered_aplakidou.faa") if l[0] == ">"]
qset = set(queries)
sub = [l.strip() for l in open(R + "/nmpfams2_10k_queries.txt")]

fam50 = collections.defaultdict(list)
for l in open(R + "/unk_clusters_50_n88.tsv"):
    a, b = l.rstrip("\n").split("\t")
    fam50[a].append(b)
fam50 = {k: v for k, v in fam50.items() if len(v) >= 3}
fam30 = collections.defaultdict(list)
for l in open(R + "/unk_clusters_30.tsv"):
    a, b = l.rstrip("\n").split("\t")
    fam30[a].append(b)
big30 = set()
for i, l in enumerate(open(R + "/unk_family_sizes_30_raw.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    if int(f[3]) >= 25: big30.add(f[0])

rows = [("source", "set", "criterion", "n_total", "n_queried", "represented", "not_represented", "pct_represented")]
def fams(src, name, crit, fam, hit):
    q = [k for k, v in fam.items() if any(m in qset for m in v)]
    r = sum(1 for k in q if any(m in hit for m in fam[k]))
    rows.append((src, name, crit, len(fam), len(q), r, len(q) - r, f"{100.0*r/len(q):.2f}" if q else "NA"))

for src in ("metag", "iso"):
    h = R + f"/nmpfams2_{src}_rep_hits.tsv"
    s = R + f"/nmpfams2_{src}_all_hits_10k.tsv"
    if not os.path.exists(h): continue
    hit = {l.split("\t", 1)[0] for l in open(h)}
    n = len(queries); r = len(hit & qset)
    rows.append((src, "proteins (95% reps, Aplakidou-filtered)", "vs family representatives", n, n, r, n - r, f"{100.0*r/n:.2f}"))
    fams(src, "50% families (unk_clusters_50_n88)", ">=3 members (95% reps)", fam50, hit)
    fams(src, "30% families (unk_clusters_30)", ">=25 complete >=35 aa proteins", {k: fam30[k] for k in big30}, hit)
    if os.path.exists(s):
        allhit = {l.split("\t", 1)[0] for l in open(s)}
        a = sum(1 for k in sub if k in allhit); b = sum(1 for k in sub if k in hit)
        extra = sum(1 for k in sub if k in allhit and k not in hit)
        lost = sum(1 for k in sub if k in hit and k not in allhit)
        rows.append((src, "10k random queries (seed 26)", "vs representatives", len(sub), len(sub), b, len(sub) - b, f"{100.0*b/len(sub):.2f}"))
        rows.append((src, "10k random queries (seed 26)", "vs all members", len(sub), len(sub), a, len(sub) - a, f"{100.0*a/len(sub):.2f}"))
        rows.append((src, "10k random queries (seed 26)", "hit vs all members only (extra)", len(sub), len(sub), extra, "NA", f"{100.0*extra/len(sub):.2f}"))
        rows.append((src, "10k random queries (seed 26)", "hit vs representatives only", len(sub), len(sub), lost, "NA", f"{100.0*lost/len(sub):.2f}"))

with open(OUT, "w") as fh:
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
for r in rows: print("\t".join(map(str, r)))
print("wrote", OUT)
