# Step 34j A0: members and strict eligibility labels of the 4,509 known length-matched controls (K reps,
#   results/step34d_known_lenmatched_fixed.faa), with the Step 34g contig rule (scripts/step34g_prep.py, unchanged):
#   a member contig takes the label of any supported rep on it (results/step34a_rep_origin_v2_final.tsv); otherwise
#   Tiara call (results/tiara/<stem>.tiara.tsv; none = < 1 kb, not eligible) and UniRef50 LCA v2 domain; unlabelled
#   Tiara-prokaryotic non-Eukaryota contigs need minimap2 -c to dm6 and the frog. Contigs already mapped in Step 34g
#   (results/step34g/maplists) reuse those PAFs; the rest are listed in results/step34j/ctl_maplists/<stem>.txt.
# Outputs (refuse to overwrite): results/step34j/ctl_members.tsv (columns as step34g/members.tsv plus paf_dir),
#   results/step34j/ctl_stems_to_map.txt, results/step34j/ctl_maplists/.
import os, gzip, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34j/"; G = R + "step34g/"
LCA = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_scaffold_classify_mmseqs"
if os.path.exists(D + "ctl_members.tsv"): raise SystemExit("REFUSING: ctl_members.tsv exists")
os.makedirs(D + "ctl_maplists")
GRP = {"bacteria": "Tiara-prokaryotic", "archaea": "Tiara-prokaryotic", "prokarya": "Tiara-prokaryotic",
       "eukarya": "Tiara-eukaryotic other", "organelle": "Tiara-eukaryotic other", "unknown": "Tiara-unknown"}
def contig(g):
    s, rest = g.split("__", 1); return (s, rest.rsplit("_", 1)[0])
def ldom(lin):
    for t in lin.split(";"):
        n = t.split("_", 1)[-1]
        if n in ("Bacteria", "Archaea", "Eukaryota", "Viruses"): return n
    return "unclassified"
K = {l[1:].split()[0] for l in open(R + "step34d_known_lenmatched_fixed.faa") if l[0] == ">"}
if len(K) != 4509: raise SystemExit("ABORT: controls %d" % len(K))
mem = []
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if r in K: mem.append((r, m))
print("control member proteins %d | reps with members %d" % (len(mem), len({r for r, _ in mem})), flush=True)
if len({r for r, _ in mem}) != 4509: raise SystemExit("ABORT: control reps missing from clusters_95_n88")
mc = {contig(m) for _, m in mem}
lab = collections.defaultdict(set)
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ix = {c: h.index(c) for c in h}
    for l in fh:
        f = l.rstrip("\n").split("\t"); k = contig(f[0])
        if k in mc: lab[k].add((f[ix["origin"]], f[ix["fly90"]], f[ix["lca_domain_v2"]]))
bad = [k for k, v in lab.items() if len(v) > 1]
if bad: raise SystemExit("ABORT: %d contigs with inconsistent rep labels" % len(bad))
need = {k for k in mc if k not in lab}
print("member contigs %d | with a supported-rep label %d | without %d" % (len(mc), len(lab), len(need)), flush=True)
stems = sorted({k[0] for k in need}); tia, dom = {}, {}
for s in stems:
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 2); c = f[0].split()[0]
            if (s, c) in need: tia[(s, c)] = GRP[f[1]]
    p = "%s/%s/%s_uniref50_lca.tsv" % (LCA, s, s)
    if not os.path.exists(p): p += ".gz"
    if not os.path.exists(p): continue
    for l in (gzip.open(p, "rt") if p.endswith(".gz") else open(p)):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 10 and (s, f[0]) in need: dom[(s, f[0])] = ldom(f[9])
done = set()
for s in (l.strip() for l in open(G + "stems_to_map.txt") if l.strip()):
    for c in open(G + "maplists/%s.txt" % s):
        if c.strip(): done.add((s, c.strip()))
tomap = collections.defaultdict(list); reuse = 0
for k in need:
    if tia.get(k) == "Tiara-prokaryotic" and dom.get(k) != "Eukaryota":
        if k in done: reuse += 1
        else: tomap[k[0]].append(k[1])
for s, v in tomap.items():
    with open(D + "ctl_maplists/%s.txt" % s, "w") as o: o.write("\n".join(sorted(v)) + "\n")
with open(D + "ctl_stems_to_map.txt", "w") as o: o.write("\n".join(sorted(tomap)) + "\n")
print("unlabelled contigs: Tiara-prokaryotic %d, other Tiara %d, no call %d | to map: reuse step34g PAF %d, new %d in %d stems"
      % (sum(1 for k in need if tia.get(k) == "Tiara-prokaryotic"), sum(1 for k in need if k in tia and tia[k] != "Tiara-prokaryotic"),
         sum(1 for k in need if k not in tia), reuse, sum(map(len, tomap.values())), len(tomap)), flush=True)
with open(D + "ctl_members.tsv", "w") as o:
    o.write("rep\tmember\tstem\tcontig\tlabel_source\torigin\tfly90\tlca_domain_v2\tpending_map\tpaf_dir\n")
    for r, m in mem:
        k = contig(m)
        if k in lab: a, b, c = list(lab[k])[0]; src = "step34a_rep_label"; pend = 0; pd = "-"
        else:
            a = tia.get(k, "<1 kb unclassified"); b = "NA"; c = dom.get(k, "NA"); src = "computed"
            pend = int(a == "Tiara-prokaryotic" and c != "Eukaryota")
            pd = ("step34g" if k in done else "step34j") if pend else "-"
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n" % (r, m, k[0], k[1], src, a, b, c, pend, pd))
print("wrote", D + "ctl_members.tsv")
