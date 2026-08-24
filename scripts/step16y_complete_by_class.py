import gzip, glob, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
cls={}
with open(W+"/results/rep_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if p[1]=="1": cls[p[0]]=p[5]
print("supported reps:",len(cls),flush=True)
L={}
for ln in open(W+"/results/supported_cds_lengths.tsv"):
    g,l=ln.rstrip("\n").split("\t"); L[g]=int(l)

OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
stems=(open(W+"/metadata/wf22_stems_44.txt").read().split()
     + open(W+"/metadata/new45_stems.txt").read().split())
print("stems:",len(stems),flush=True)
code={}
found=0
for i,s in enumerate(stems,1):
    p=NEW+"/"+s+".aa.fa.gz"
    try: fh=gzip.open(p,"rt")
    except: fh=gzip.open(OLD+"/"+s+".aa.fa.gz","rt")
    for line in fh:
        if line[0]!=">": continue
        tok=line[1:].split("#")
        gid=s+"__"+tok[0].split()[0]
        if gid not in cls: continue
        j=line.find("partial=")
        if j>=0:
            code[gid]=line[j+8:j+10]; found+=1
    fh.close()
    if i%20==0: print("stems done:",i,"found:",found,flush=True)
print("reps with a flag:",len(code),"of",len(cls))

lab={"00":"complete","10":"no start","01":"no stop","11":"both missing"}
tab=collections.defaultdict(collections.Counter)
lensum=collections.defaultdict(int); lenn=collections.defaultdict(int)
for g,k in cls.items():
    c=code.get(g,"NA")
    tab[k][c]+=1
    if c=="00":
        lensum[k]+=L.get(g,0); lenn[k]+=1
print("\n=== ORF completeness by class ===")
for k in ["K","KWP","U"]:
    t=sum(tab[k].values())
    parts=" ".join("%s:%d(%.1f%%)"%(lab.get(c,c),tab[k][c],100.0*tab[k][c]/t)
                   for c in ["00","10","01","11","NA"] if tab[k][c])
    print(k,"n=%d"%t,parts)
print("\n=== mean length of COMPLETE genes only ===")
for k in ["K","KWP","U"]:
    if lenn[k]: print("%s complete n=%d mean_bp=%.0f"%(k,lenn[k],lensum[k]/lenn[k]))
with open(W+"/results/completeness_by_class.tsv","w") as fh:
    fh.write("class\tcode\tlabel\tn\n")
    for k in ["K","KWP","U"]:
        for c in tab[k]: fh.write("%s\t%s\t%s\t%d\n"%(k,c,lab.get(c,c),tab[k][c]))
print("\nwrote results/completeness_by_class.tsv")
