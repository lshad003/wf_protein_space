# Step 24b: gene accumulation vs number of samples, from the count matrix.
# 100 permutations of sample order; reports mean and SD of cumulative
# genes detected, split known (K+KWP) vs unknown (GU+EU), plus all.
import gzip, random, numpy as np
W="/bigdata/stajichlab/lshad003/wf_protein_space"
random.seed(1); np.random.seed(1)

cls={}
with open(W+"/results/fourway_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t"); cls[p[0]]=p[1]
print("classes loaded:",len(cls),flush=True)

fh=gzip.open(W+"/results/count_matrix_primary.tsv.gz","rt")
hdr=fh.readline().rstrip("\n").split("\t")[1:]
ns=len(hdr)
grp={"all":[], "known":[], "unknown":[]}
rows=[]
n=0
for ln in fh:
    p=ln.rstrip("\n").split("\t")
    v=np.array([1 if int(x)>0 else 0 for x in p[1:]],dtype=np.uint8)
    if v.sum()==0: continue
    c=cls.get(p[0])
    rows.append((v, "known" if c in ("K","KWP") else "unknown"))
    n+=1
    if n%400000==0: print("rows:",n,flush=True)
print("genes with any count:",len(rows),flush=True)

M=np.array([r[0] for r in rows],dtype=np.uint8)
lab=np.array([r[1] for r in rows])
del rows
idx={"all":np.ones(len(lab),bool),"known":lab=="known","unknown":lab=="unknown"}

PERM=100
out={k:np.zeros((PERM,ns)) for k in idx}
order=list(range(ns))
for p in range(PERM):
    random.shuffle(order)
    for k,m in idx.items():
        sub=M[m]
        seen=np.zeros(sub.shape[0],bool)
        for i,s in enumerate(order):
            seen |= sub[:,s].astype(bool)
            out[k][p,i]=seen.sum()
    if (p+1)%20==0: print("perm",p+1,flush=True)

with open(W+"/results/accumulation_by_sample.tsv","w") as f:
    f.write("set\tn_samples\tmean\tsd\n")
    for k in ["all","known","unknown"]:
        for i in range(ns):
            f.write("%s\t%d\t%.1f\t%.1f\n"%(k,i+1,out[k][:,i].mean(),out[k][:,i].std()))
print("\nwrote results/accumulation_by_sample.tsv")
for k in ["all","known","unknown"]:
    a=out[k].mean(axis=0)
    print("%-8s 1 sample %.0f | 45 %.0f | 89 %.0f | last 10 samples add %.2f%%"%(
        k,a[0],a[44],a[88],100*(a[88]-a[78])/a[88]))
