# Step 20c: taxonomy of each supported gene from its contig's LCA assignment.
# Gene IDs are stem__contig_gene, so the contig is recoverable from the ID.
# Coverage limited to the stems that have a uniref50_lca.tsv file.
import os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
FEC="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_scaffold_classify_mmseqs"

cls={}
with open(W+"/results/fourway_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t"); cls[p[0]]=p[1]
print("supported genes:",len(cls),flush=True)

stems=(open(W+"/metadata/wf22_stems_44.txt").read().split()
     + open(W+"/metadata/wf23_stems.txt").read().split()
     + open(W+"/metadata/wf24_stems.txt").read().split())

tax={}
have=0
for i,s in enumerate(stems,1):
    f="%s/%s/%s_uniref50_lca.tsv"%(FEC,s,s)
    if not os.path.exists(f): continue
    have+=1
    for ln in open(f):
        p=ln.rstrip("\n").split("\t")
        if len(p)<10: continue
        tax[s+"__"+p[0]] = (p[2], p[3], p[9])
    if i%20==0: print("stems",i,"with tax",have,flush=True)
print("stems with taxonomy:",have,"of",len(stems))
print("contigs with taxonomy:",len(tax),flush=True)

def contig_of(g):
    stem,rest = g.split("__",1)
    parts = rest.rsplit("_",1)
    return stem+"__"+parts[0]

def domain(lin):
    for t in lin.split(";"):
        if t.startswith("d_"): return t[2:]
    return "unclassified"

out=open(W+"/results/gene_taxonomy.tsv","w")
out.write("rep\tclass4\trank\tname\tdomain\n")
dom=collections.defaultdict(collections.Counter)
notax=collections.Counter()
for g,c in cls.items():
    t=tax.get(contig_of(g))
    if t is None:
        notax[c]+=1; continue
    d=domain(t[2])
    dom[c][d]+=1
    out.write("%s\t%s\t%s\t%s\t%s\n"%(g,c,t[0],t[1],d))
out.close()

print("\n=== domain assignment by class (of genes with a classified contig) ===")
for c in ["K","KWP","GU","EU"]:
    tot=sum(dom[c].values())
    if not tot: continue
    top=" ".join("%s:%.1f%%"%(k,100.0*v/tot) for k,v in dom[c].most_common(5))
    print("%-5s n=%-9d %s"%(c,tot,top))
print("\ngenes with no contig taxonomy:", dict(notax))
print("wrote results/gene_taxonomy.tsv")
