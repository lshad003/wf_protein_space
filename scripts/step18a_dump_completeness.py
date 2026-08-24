import gzip
W="/bigdata/stajichlab/lshad003/wf_protein_space"
keep=set()
with open(W+"/results/rep_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if p[1]=="1": keep.add(p[0])
OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
stems=(open(W+"/metadata/wf22_stems_44.txt").read().split()
     + open(W+"/metadata/new45_stems.txt").read().split())
out=open(W+"/results/completeness_reps.tsv","w"); n=0
for i,s in enumerate(stems,1):
    try: fh=gzip.open(NEW+"/"+s+".aa.fa.gz","rt")
    except: fh=gzip.open(OLD+"/"+s+".aa.fa.gz","rt")
    for line in fh:
        if line[0]!=">": continue
        gid=s+"__"+line[1:].split("#")[0].split()[0]
        if gid not in keep: continue
        j=line.find("partial=")
        if j>=0:
            out.write("%s\t%s\n"%(gid,line[j+8:j+10])); n+=1
    fh.close()
    if i%20==0: print("stems",i,"rows",n,flush=True)
out.close()
print("wrote",n,"expect 2069453")
