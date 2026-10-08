# Step 16w n88: class abundance shares on results/count_matrix_primary_n88.tsv.gz
# with classes from results/rep_classes_n88.tsv.
import gzip, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
cls={}
with open(W+"/results/rep_classes_n88.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if p[1]=="1": cls[p[0]]=p[5]
L={}
for ln in open(W+"/results/supported_cds_lengths.tsv"):
    g,l=ln.rstrip("\n").split("\t"); L[g]=int(l)
coh={}
for ln in open(W+"/results/cohort_map.tsv"):
    s,c=ln.rstrip("\n").split("\t"); coh[s]=c
fh=gzip.open(W+"/results/count_matrix_primary_n88.tsv.gz","rt")
hdr=fh.readline().rstrip("\n").split("\t")[1:]
ns=len(hdr)
raw={k:[0]*ns for k in ["K","KWP","U"]}
rpk={k:[0.0]*ns for k in ["K","KWP","U"]}
n=0
for ln in fh:
    p=ln.rstrip("\n").split("\t")
    k=cls.get(p[0]); 
    if k is None: continue
    ln_kb=L[p[0]]/1000.0
    r=raw[k]; q=rpk[k]
    for j in range(ns):
        x=int(p[j+1])
        if x:
            r[j]+=x
            q[j]+=x/ln_kb
    n+=1
    if n%500000==0: print(n,flush=True)
out=open(W+"/results/dark_abundance_by_sample_n88.tsv","w")
out.write("stem\tcohort\tpct_reads_K\tpct_reads_KWP\tpct_reads_U\tpct_rpk_U\n")
agg=collections.defaultdict(list)
for j,s in enumerate(hdr):
    t=raw["K"][j]+raw["KWP"][j]+raw["U"][j]
    tq=rpk["K"][j]+rpk["KWP"][j]+rpk["U"][j]
    vals=(100.0*raw["K"][j]/t,100.0*raw["KWP"][j]/t,100.0*raw["U"][j]/t,100.0*rpk["U"][j]/tq)
    out.write("%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n"%((s,coh[s])+vals))
    agg[coh[s]].append(vals[2])
out.close()
print("=== percent of mapped reads going to unknown (U) genes ===")
for c in sorted(agg):
    v=sorted(agg[c]); m=sum(v)/len(v)
    print("%s n=%d mean=%.2f%% min=%.2f%% max=%.2f%%"%(c,len(v),m,v[0],v[-1]))
allv=[x for c in agg for x in agg[c]]
print("all %d: mean=%.2f%%  per-sample range %.2f%% to %.2f%%"%(len(allv),sum(allv)/len(allv),min(allv),max(allv)))
nU=sum(1 for k in cls.values() if k=="U")
print("pooled raw reads: U %.2f%% of mapped reads to classed genes; length-normalized U %.2f%%"%(100.0*sum(raw["U"])/sum(sum(raw[k]) for k in raw), 100.0*sum(rpk["U"])/sum(sum(rpk[k]) for k in rpk)))
print("reference: U is %.2f%% of supported genes (%d of %d)"%(100.0*nU/len(cls),nU,len(cls)))
