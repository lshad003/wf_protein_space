# Step 34o: the strict_final representatives present in no cohort under the Fig. 1c presence rule (Step 34k,
#   scripts/step34k_final.py:6-8,98-106: present in a sample = a 95% cluster member from that sample on a strict-eligible
#   contig, results/step34k/contig_labels.tsv strict_eligible == 1; the 88 samples of metadata/samples_88_treatment.tsv).
# Expected 91 (results/step34k_summary.tsv, "strict_members none (all classes)").
# Per representative: class4 and origin (results/step34a_rep_origin_v2_final.tsv), source sample and contig, members in
#   the 88 samples (results/clusters_95_n88.tsv), and for each member contig its Step 34k label and length (len= in the
#   MEGAHIT headers of <FEC>/<stem>/scaffolds/<stem>_R.fa.gz, as step32a_tiara_array.sh).
# Output: results/step34o_nocohort_reps.tsv; refuses to overwrite.
import os, gzip, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
FEC = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results"
OUT = R + "step34o_nocohort_reps.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
stems = {l.split("\t")[0] for l in list(open(W + "/metadata/samples_88_treatment.tsv"))[1:]}
if len(stems) != 88: raise SystemExit("ABORT: %d samples" % len(stems))
cls = {}; ori = {}
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ic, isf, io = h.index("class4"), h.index("strict_final"), h.index("origin")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[isf] == "1": cls[f[0]] = f[ic]; ori[f[0]] = f[io]
if len(cls) != 1359273: raise SystemExit("ABORT: %d strict_final reps" % len(cls))
lab = {}
with open(R + "step34k/contig_labels.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t"))); lab[f["stem"] + "\t" + f["contig"]] = f
def ckey(m):
    s, rest = m.split("__", 1); return s, rest.rsplit("_", 1)[0]
elig = set(); seen = set()
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if r not in cls: continue
        s, c = ckey(m)
        if s not in stems: continue
        seen.add(r)
        x = lab.get(s + "\t" + c)
        if x is not None and x["strict_eligible"] == "1": elig.add(r)
none = sorted(r for r in seen if r not in elig)
print("strict_final reps with a member in the 88: %d; with no eligible member: %d" % (len(seen), len(none)))
if len(none) != 91: raise SystemExit("ABORT: expected 91")
mem = collections.defaultdict(list); want = collections.defaultdict(set); ns = set(none)
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if r not in ns: continue
        s, c = ckey(m)
        if s in stems: mem[r].append((s, c)); want[s].add(c)
for r in none:
    s, c = ckey(r); want[s].add(c)
clen = {}
for s, cs in sorted(want.items()):
    p = "%s/%s/scaffolds/%s_R.fa.gz" % (FEC, s, s)
    with gzip.open(p, "rt") as fh:
        for l in fh:
            if l[0] != ">": continue
            f = l[1:].split(); 
            if f[0] in cs: clen[(s, f[0])] = int([x for x in f if x.startswith("len=")][0][4:])
miss = sum(1 for s, cs in want.items() for c in cs if (s, c) not in clen)
print("contig lengths found %d, missing %d" % (len(clen), miss))
with open(OUT, "w") as o:
    o.write("rep\tclass4\trep_origin\tsource_sample\trep_contig\trep_contig_len\tn_members_in_88\tn_member_samples\t"
            "max_member_contig_len\tmember_contig_origins\tmember_contigs\n")
    for r in none:
        s, c = ckey(r); M = mem[r]
        ls = [clen.get(x, -1) for x in M]
        oc = collections.Counter(lab[a + "\t" + b]["origin"] if (a + "\t" + b) in lab else "unlabelled" for a, b in M)
        det = ";".join("%s|%s|%s|%s|%s" % (a, b, clen.get((a, b), "NA"), lab.get(a + "\t" + b, {}).get("origin", "unlabelled"),
                       lab.get(a + "\t" + b, {}).get("strict_eligible", "NA")) for a, b in M)
        o.write("\t".join(map(str, (r, cls[r], ori[r], s, c, clen.get((s, c), "NA"), len(M), len({a for a, b in M}),
                max(ls), ",".join("%s:%d" % kv for kv in oc.most_common()), det))) + "\n")
print("by class:", dict(collections.Counter(cls[r] for r in none)))
print("source sample:", dict(collections.Counter(ckey(r)[0] for r in none)))
print("members: total %d, max member contig length %d" % (sum(len(mem[r]) for r in none), max(clen.get(x, -1) for r in none for x in mem[r])))
print("wrote", OUT)
