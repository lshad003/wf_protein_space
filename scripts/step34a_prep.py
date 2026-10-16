# Step 34a prep: contigs that carry a supported 95% rep (results/fourway_classes_n88.tsv, 2,067,011 reps, n88 catalog
# incl. reps on UHM586.41010 contigs) and were not mapped in step33c2 / step33g (i.e. not in results/step33c_lists).
# Protein -> contig parsed as step29d. Outputs (refuse to overwrite): results/step34a_lists/<stem>.txt (one contig per
# line) and results/step34a_lists/stems.txt (array order: the 88 of metadata/samples_88_treatment.tsv, then UHM586.41010).
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; LD = R + "step34a_lists"
if os.path.exists(LD): raise SystemExit("REFUSING: " + LD + " exists")
X = "UHM586.41010"
stems = [l.split("\t")[0] for i, l in enumerate(open(W + "/metadata/samples_88_treatment.tsv")) if i] + [X]
done = set()
for s in stems:
    for l in open(R + "step33c_lists/%s.tsv" % s): done.add((s, l.split("\t")[0]))
rc = set(); nrep = 0
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        r = l.split("\t", 1)[0]; s, rest = r.split("__", 1); nrep += 1
        if s not in stems: raise SystemExit("ABORT: rep from stem " + s)
        rc.add((s, rest[:rest.rfind("_")]))
if nrep != 2067011: raise SystemExit("ABORT: %d supported reps" % nrep)
todo = rc - done; per = collections.defaultdict(list)
for s, c in todo: per[s].append(c)
os.mkdir(LD)
for s in stems:
    with open(LD + "/%s.txt" % s, "w") as o:
        for c in sorted(per[s]): o.write(c + "\n")
open(LD + "/stems.txt", "w").write("\n".join(stems) + "\n")
print("supported reps %d on %d contigs | already mapped (step33c lists) %d | to map now %d" % (nrep, len(rc), len(rc & done), len(todo)))
print("contigs to map per stem: min %d, median %d, max %d" % (min(len(per[s]) for s in stems), sorted(len(per[s]) for s in stems)[len(stems) // 2], max(len(per[s]) for s in stems)))
