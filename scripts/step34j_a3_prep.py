# Step 34j A3 prep: contigs and proteins for geNomad and DefenseFinder.
# Contigs = every contig carrying a strict member of a recurrent EU family (results/step34j/context_members.tsv,
#   set EU; the 186 structural-hit families are a subset and are reported separately downstream).
# Contig sequences from the assemblies (shared Fecal/results/<stem>/scaffolds/<stem>_R.fa.gz), renamed <stem>__<contig>.
# Proteins of every gene on those contigs from the per-assembly Prodigal files (catalog/by_assembly_new, else shared
#   Proteins/by_assembly), named <stem>__<contig>_<n> (catalog IDs), terminal '*' stripped, written contig by contig
#   in gene order (gembase layout: replicon = ID before the last '_'), for DefenseFinder --db-type gembase.
# Outputs (refuse to overwrite): results/step34j/mge/contigs.fna, results/step34j/mge/proteins.faa
import os, gzip, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; D = W + "/results/step34j/"; O = D + "mge/"
FEC = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results"
OLD = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"; NEW = W + "/catalog/by_assembly_new"
if os.path.exists(O): raise SystemExit("REFUSING: " + O + " exists")
tc = collections.defaultdict(set)
with open(D + "context_members.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); i_s, i_st, i_c = h.index("set"), h.index("stem"), h.index("contig")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[i_s] == "EU": tc[f[i_st]].add(f[i_c])
print("contigs %d in %d stems" % (sum(map(len, tc.values())), len(tc)), flush=True)
os.makedirs(O); nc = np_ = 0
with open(O + "contigs.fna", "w") as oc, open(O + "proteins.faa", "w") as op:
    for s in sorted(tc):
        want = tc[s]; w = False; got = set()
        with gzip.open("%s/%s/scaffolds/%s_R.fa.gz" % (FEC, s, s), "rt") as fh:
            for l in fh:
                if l[0] == ">":
                    c = l[1:].split()[0]; w = c in want
                    if w: oc.write(">%s__%s\n" % (s, c)); got.add(c); nc += 1
                elif w: oc.write(l)
        if got != want: raise SystemExit("ABORT: %s: %d contigs not found in the assembly" % (s, len(want - got)))
        p = NEW + "/" + s + ".aa.fa.gz"
        if not os.path.exists(p): p = OLD + "/" + s + ".aa.fa.gz"
        prot = collections.defaultdict(dict); cur = None
        with gzip.open(p, "rt") as fh:
            for l in fh:
                if l[0] == ">":
                    gid = l[1:].split(" # ")[0].strip(); c, n = gid.rsplit("_", 1)
                    cur = (c, int(n)) if c in want else None
                    if cur: prot[c][int(n)] = []
                elif cur: prot[cur[0]][cur[1]].append(l.strip())
        for c in sorted(prot):
            for n in sorted(prot[c]):
                op.write(">%s__%s_%d\n%s\n" % (s, c, n, "".join(prot[c][n]).rstrip("*"))); np_ += 1
print("wrote contigs %d, proteins %d" % (nc, np_))
