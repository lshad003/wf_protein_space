import gzip, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
W="/bigdata/stajichlab/lshad003/wf_protein_space"
C={"K":"#0072B2","KWP":"#E69F00","U":"#D55E00",
   "WF22":"#0072B2","WF23":"#009E73","WF24":"#CC79A7"}
BINS=["1","2-5","6-20","21-60","61-88","89"]
def b(x):
    if x==1: return "1"
    if x<=5: return "2-5"
    if x<=20: return "6-20"
    if x<=60: return "21-60"
    if x<89: return "61-88"
    return "89"
agg=collections.defaultdict(collections.Counter)
with gzip.open(W+"/results/prevalence_primary.tsv.gz","rt") as fh:
    next(fh)
    for ln in fh:
        p=ln.split("\t")
        agg[p[1]][b(int(p[2]))]+=1
fig,ax=plt.subplots(figsize=(7,4.2))
w=0.26
for i,k in enumerate(["K","KWP","U"]):
    t=sum(agg[k].values())
    y=[100.0*agg[k][x]/t for x in BINS]
    ax.bar([j+(i-1)*w for j in range(len(BINS))],y,w,label=k,color=C[k])
ax.set_xticks(range(len(BINS))); ax.set_xticklabels(BINS)
ax.set_xlabel("Number of samples the gene is detected in (of 89)")
ax.set_ylabel("Percent of genes in class")
ax.set_title("Prevalence of supported genes by annotation class")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(W+"/results/fig_prevalence_by_class.png",dpi=300)
print("fig1 written")

rows=[l.rstrip("\n").split("\t") for l in open(W+"/results/dark_abundance_by_sample.tsv")][1:]
fig,ax=plt.subplots(figsize=(7,4.2))
xs={"WF22":0,"WF23":1,"WF24":2}
import random; random.seed(1)
for r in rows:
    x=xs[r[1]]; j=(random.random()-0.5)*0.28
    ax.scatter(x-0.18+j,float(r[4]),color=C[r[1]],alpha=0.75,s=22,marker="o")
    ax.scatter(x+0.18+j,float(r[5]),color=C[r[1]],alpha=0.75,s=22,marker="^")
ax.set_xticks([0,1,2]); ax.set_xticklabels(["WF22 (n=44)","WF23 (n=9)","WF24 (n=36)"])
ax.set_ylabel("Percent of gut gene abundance")
ax.set_title("Unannotated (U) share per sample: circles raw reads, triangles length-normalized")
ax.axhline(7.03,ls=":",c="grey",lw=1); ax.axhline(19.28,ls="--",c="grey",lw=1)
fig.tight_layout(); fig.savefig(W+"/results/fig_dark_abundance.png",dpi=300)
print("fig2 written")
