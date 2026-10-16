# Step 33c prep: per-sample contig lists for host (wood frog) mapping, and the combined Basidiobolus reference.
# Sets (contig keys = stem + first word of the contig header; protein -> contig parsed as step29d):
#  A  Tiara-eukaryotic contigs >= 1 kb (eukarya + organelle; results/tiara, 88 samples)
#  B  contigs carrying an EU supported rep (results/fourway_classes_n88.tsv class4 EU; any length). As in step32b,
#     the 405 EU reps whose own contig is in the excluded control UHM586.41010 are included (B only).
#  C  random 2% of Tiara-prokaryotic contigs >= 3 kb (bacteria + archaea + prokarya; 88 samples),
#     random.Random(33).sample over the sorted contig keys; negative control
#  F  contigs carrying a member of the 15,580 recurrent EU families (results/unk_direct30_clusters.tsv, members from
#     the 88 samples only); needed for the per-family host majority
# Outputs (refuse to overwrite): results/step33c_lists/<stem>.tsv (contig, sets), results/step33c_lists/stems.txt
#  (array order: the 88 of metadata/samples_88_treatment.tsv, then UHM586.41010),
#  catalog/basid/basid6.fa (six genomes of FINAL_FROZEN_GENOMES_BASIDIOBOLUS_AUG18_2025, headers "<strain>|<scaffold>")
import os, re, gzip, random, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
LD = R + "step33c_lists"; BF = W + "/catalog/basid/basid6.fa"
BD = "/bigdata/stajichlab/shared/projects/Herptile/Basidiobolus/FINAL_FROZEN_GENOMES_BASIDIOBOLUS_AUG18_2025"
STRAINS = (("STP1710", "Basidiobolus_STP1710_7"), ("STP1717", "Basidiobolus_STP1717_1"), ("UHM207", "Basidiobolus_UHM207_4505"),
           ("UHM260", "Basidiobolus_UHM260_5136"), ("UHM516", "Basidiobolus_UHM516_7697"), ("UHM520", "Basidiobolus_UHM520_7734"))
X = "UHM586.41010"
for p in (LD, BF):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
stems = [l.split("\t")[0] for i, l in enumerate(open(W + "/metadata/samples_88_treatment.tsv")) if i]
if len(stems) != 88: raise SystemExit("ABORT: %d stems" % len(stems))
S88 = set(stems)
EUK, PROK = {"eukarya", "organelle"}, {"bacteria", "archaea", "prokarya"}
sets = collections.defaultdict(set); prok3 = []
for s in stems:
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 2); c = f[0].split()[0]; ln = int(re.search(r"len=(\d+)", f[0]).group(1))
            if f[1] in EUK: sets[(s, c)].add("A")
            elif f[1] in PROK and ln >= 3000: prok3.append((s, c))
nA = sum(1 for v in sets.values() if "A" in v)
prok3.sort(); kC = round(0.02 * len(prok3))
for k in random.Random(33).sample(prok3, kC): sets[k].add("C")
def contig(g):
    s, rest = g.split("__", 1); return s, rest[:rest.rfind("_")]
nB = collections.Counter()
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        r, c = l.split("\t", 2)[:2]
        if c != "EU": continue
        k = contig(r)
        if k[0] not in S88 and k[0] != X: raise SystemExit("ABORT: EU rep from stem " + k[0])
        sets[k].add("B"); nB[k[0] == X] += 1
fam = set()
with open(R + "unk_direct30_eu_recurrence.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ia, ic = h.index("n_animals"), h.index("n_cohorts")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if int(f[ia]) >= 10 and int(f[ic]) >= 2: fam.add(f[0])
if len(fam) != 15580: raise SystemExit("ABORT: %d families" % len(fam))
nF = 0
for l in open(R + "unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in fam:
        k = contig(m)
        if k[0] not in S88: raise SystemExit("ABORT: family member from stem " + k[0])
        sets[k].add("F"); nF += 1
os.mkdir(LD); per = collections.defaultdict(list)
for (s, c), v in sets.items(): per[s].append((c, "".join(sorted(v))))
order = stems + [X]
for s in order:
    with open(LD + "/%s.tsv" % s, "w") as o:
        for c, v in sorted(per[s]): o.write("%s\t%s\n" % (c, v))
open(LD + "/stems.txt", "w").write("\n".join(order) + "\n")
cnt = collections.Counter(x for v in sets.values() for x in v)
print("A Tiara-eukaryotic >= 1 kb contigs: %d" % nA)
print("C random 2%% of %d Tiara-prokaryotic >= 3 kb contigs (seed 33): %d" % (len(prok3), kC))
print("B EU supported reps: %d (of them on %s: %d) -> contigs %d" % (sum(nB.values()), X, nB[True], cnt["B"]))
print("F family members: %d -> contigs %d" % (nF, cnt["F"]))
print("contigs in any set: %d | overlaps A&B %d, B&F %d, A&F %d, C&B %d" % (len(sets), sum(1 for v in sets.values() if {"A", "B"} <= v),
      sum(1 for v in sets.values() if {"B", "F"} <= v), sum(1 for v in sets.values() if {"A", "F"} <= v), sum(1 for v in sets.values() if {"B", "C"} <= v)))
print("contigs per task: min %d, max %d (%s: %d)" % (min(len(per[s]) for s in order), max(len(per[s]) for s in order), X, len(per[X])))
os.makedirs(os.path.dirname(BF), exist_ok=True); n = 0
with open(BF, "w") as o:
    for st, pre in STRAINS:
        for l in gzip.open("%s/%s_UHMproject_genome_FINAL/%s_softmasked.fasta.gz" % (BD, st, pre), "rt"):
            if l[0] == ">": o.write(">%s|%s\n" % (st, l[1:].split()[0])); n += 1
            else: o.write(l)
print("wrote %s: %d scaffolds from 6 strains" % (BF, n))
