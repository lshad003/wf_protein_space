# Figure: unannotated genes are small proteins, not fragments.
# A length distributions by class, complete ORFs only.
# B mean length all vs complete only, the key comparison.
# C AntiFam rate. D length by number of years detected.
import collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
W="/bigdata/stajichlab/lshad003/wf_protein_space"
C={"K":"#0072B2","KWP":"#56B4E9","GU":"#E69F00","EU":"#D55E00"}
O=["K","KWP","GU","EU"]

L=collections.defaultdict(list); Lc=collections.defaultdict(list)
sp=collections.Counter(); n=collections.Counter()
yl=collections.defaultdict(list)
with open(W+"/results/fourway_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        c=p[1]; b=int(p[2]); n[c]+=1
        L[c].append(b)
        if p[3]=="00": Lc[c].append(b)
        if p[7]=="1": sp[c]+=1
        if c in ("GU","EU"): yl[p[5]].append(b)

fig,ax=plt.subplots(2,2,figsize=(9,7))
for c in O:
    v=np.array(Lc[c]); v=v[v<=3000]
    ax[0,0].hist(v,bins=60,density=True,histtype="step",lw=1.6,color=C[c],label=c)
ax[0,0].set_xlabel("CDS length (bp), complete ORFs only")
ax[0,0].set_ylabel("Density"); ax[0,0].legend(frameon=False,fontsize=8)
ax[0,0].set_title("A  Length distribution",fontsize=10)

w=0.38
a=[np.mean(L[c]) for c in O]; b=[np.mean(Lc[c]) for c in O]
ax[0,1].bar([i-w/2 for i in range(4)],a,w,color="#bbbbbb",label="all genes")
ax[0,1].bar([i+w/2 for i in range(4)],b,w,color=[C[c] for c in O],label="complete ORFs only")
ax[0,1].set_xticks(range(4)); ax[0,1].set_xticklabels(O)
ax[0,1].set_ylabel("Mean CDS length (bp)"); ax[0,1].legend(frameon=False,fontsize=8)
ax[0,1].set_title("B  Excluding incomplete frames does not\nshorten the unannotated classes",fontsize=10)

v=[100.0*sp[c]/n[c] for c in O]
ax[1,0].bar(range(4),v,color=[C[c] for c in O])
ax[1,0].set_xticks(range(4)); ax[1,0].set_xticklabels(O)
ax[1,0].set_ylabel("Percent flagged by AntiFam")
ax[1,0].set_title("C  Spurious ORF screen",fontsize=10)
for i,x in enumerate(v): ax[1,0].text(i,x,"%.3f"%x,ha="center",va="bottom",fontsize=8)

ks=["3","2","1"]
m=[np.mean(yl[k]) if yl[k] else 0 for k in ks]
cnt=[len(yl[k]) for k in ks]
ax[1,1].bar(range(3),m,color="#D55E00")
ax[1,1].set_xticks(range(3))
ax[1,1].set_xticklabels(["3 years\n(n=%s)"%format(cnt[0],","),
                         "2 years\n(n=%s)"%format(cnt[1],","),
                         "1 year\n(n=%s)"%format(cnt[2],",")],fontsize=8)
ax[1,1].set_ylabel("Mean CDS length (bp)")
ax[1,1].set_title("D  Unannotated genes, by years detected",fontsize=10)
for i,x in enumerate(m): ax[1,1].text(i,x,"%.0f"%x,ha="center",va="bottom",fontsize=8)

for a2 in ax.flat:
    for s in ("top","right"): a2.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(W+"/results/fig_length_completeness.png",dpi=300)
print("wrote results/fig_length_completeness.png")
for c in O:
    print("%-4s all %.0f  complete %.0f  antifam %.3f%%"%(
        c,np.mean(L[c]),np.mean(Lc[c]),100.0*sp[c]/n[c]))
