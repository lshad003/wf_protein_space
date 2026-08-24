import gzip
W="/bigdata/stajichlab/lshad003/wf_protein_space"
cls={}
with open(W+"/results/rep_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if p[1]=="1": cls[p[0]]=p[5]
coh={}
for ln in open(W+"/results/cohort_map.tsv"):
    s,c=ln.rstrip("\n").split("\t"); coh[s]=c
fh=gzip.open(W+"/results/count_matrix_primary.tsv.gz","rt")
hdr=fh.readline().rstrip("\n").split("\t")[1:]
cvec=[coh[s] for s in hdr]
out=gzip.open(W+"/results/prevalence_primary.tsv.gz","wt")
out.write("gene\tclass\tn_samples\tn_wf22\tn_wf23\tn_wf24\ttotal\n")
import collections
dist=collections.Counter(); byc=collections.defaultdict(collections.Counter)
core=collections.Counter(); n=0
for ln in fh:
    p=ln.rstrip("\n").split("\t")
    g=p[0]; v=[int(x) for x in p[1:]]
    k=cls.get(g,"NA")
    a=b=c=0; t=0; ns=0
    for j,x in enumerate(v):
        if x>0:
            ns+=1
            if cvec[j]=="WF22": a+=1
            elif cvec[j]=="WF23": b+=1
            else: c+=1
        t+=x
    out.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n"%(g,k,ns,a,b,c,t))
    dist[ns]+=1; byc[k][ns]+=1
    if ns==89: core[k]+=1
    n+=1
    if n%500000==0: print(n,flush=True)
out.close()
print("genes:",n)
print("=== in all 89 samples, by class ===")
for k in ["K","KWP","U"]: print(k, core[k])
print("=== prevalence bins, by class ===")
def binof(x):
    if x==1: return "1"
    if x<=5: return "2-5"
    if x<=20: return "6-20"
    if x<=60: return "21-60"
    if x<89: return "61-88"
    return "89"
agg=collections.defaultdict(collections.Counter)
for k in byc:
    for ns,cnt in byc[k].items(): agg[k][binof(ns)]+=cnt
for k in ["K","KWP","U"]:
    tot=sum(agg[k].values())
    row=" ".join("%s:%d(%.1f%%)"%(b,agg[k][b],100.0*agg[k][b]/tot) for b in ["1","2-5","6-20","21-60","61-88","89"])
    print(k, "total", tot, row)
