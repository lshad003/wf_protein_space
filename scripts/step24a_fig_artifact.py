# Figure: the unannotated fraction is real, not artifact.
# a  length ECDF by class, complete ORFs only
# b  ORF completeness, stacked proportion by class
# c  spurious ORF rate by class
# d  length ECDF of unannotated genes by number of years detected
import sys, collections
sys.path.insert(0,"/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
from wfstyle import *
import numpy as np
W="/bigdata/stajichlab/lshad003/wf_protein_space"

L=collections.defaultdict(list); Lc=collections.defaultdict(list)
comp=collections.defaultdict(collections.Counter)
sp=collections.Counter(); n=collections.Counter()
yl=collections.defaultdict(list)
with open(W+"/results/fourway_classes.tsv") as fh:
    next(fh)
    for line in fh:
        p=line.rstrip("\n").split("\t")
        c=p[1]; b=int(p[2]); n[c]+=1
        L[c].append(b); comp[c][p[3]]+=1
        if p[3]=="00": Lc[c].append(b)
        if p[7]=="1": sp[c]+=1
        if c in ("GU","EU"): yl[p[5]].append(b)

def ecdf(a,x):
    a=np.sort(np.asarray(a)); return np.searchsorted(a,x,"right")/a.size

fig,ax=plt.subplots(2,2,figsize=size(180,140))
x=np.linspace(0,2500,400)
for c in CLS:
    ax[0,0].plot(x,ecdf(Lc[c],x),color=COL[c],lw=1.4,
                 label="%s (n=%s)"%(c,format(len(Lc[c]),",")))
ax[0,0].set_xlabel("CDS length (bp), complete ORFs only")
ax[0,0].set_ylabel("Cumulative fraction of genes")
ax[0,0].legend(frameon=False,loc="lower right")
ax[0,0].set_title("a   Length distribution",loc="left",fontweight="bold")

order=["00","10","01","11"]
names={"00":"complete","10":"no start","01":"no stop","11":"neither"}
shade=["#333333","#777777","#aaaaaa","#dddddd"]
bot=np.zeros(4)
for k,sh in zip(order,shade):
    v=np.array([100.0*comp[c][k]/n[c] for c in CLS])
    ax[0,1].bar(range(4),v,0.62,bottom=bot,color=sh,label=names[k])
    bot+=v
ax[0,1].set_xticks(range(4)); ax[0,1].set_xticklabels(CLS)
ax[0,1].set_ylabel("Percent of genes"); ax[0,1].set_ylim(0,100)
ax[0,1].legend(frameon=False,loc="center left",bbox_to_anchor=(1.01,0.5))
ax[0,1].set_title("b   Reading frame completeness",loc="left",fontweight="bold")

v=[100.0*sp[c]/n[c] for c in CLS]
ax[1,0].bar(range(4),v,0.62,color=[COL[c] for c in CLS])
ax[1,0].set_xticks(range(4)); ax[1,0].set_xticklabels(CLS)
ax[1,0].set_ylabel("Percent flagged, AntiFam v6.0")
ax[1,0].set_title("c   Spurious ORF screen",loc="left",fontweight="bold")
for i,q in enumerate(v):
    ax[1,0].text(i,q,"%.3f"%q if q else "0",ha="center",va="bottom",fontsize=7)

for k,lab,sh in zip(["3","2","1"],
        ["3 years","2 years","1 year"],["#D55E00","#E8916A","#F3C4AC"]):
    if yl[k]:
        ax[1,1].plot(x,ecdf(yl[k],x),color=sh,lw=1.4,
                     label="%s (n=%s)"%(lab,format(len(yl[k]),",")))
ax[1,1].set_xlabel("CDS length (bp)")
ax[1,1].set_ylabel("Cumulative fraction of genes")
ax[1,1].legend(frameon=False,loc="lower right")
ax[1,1].set_title("d   Unannotated genes, by years detected",loc="left",fontweight="bold")

fig.tight_layout()
fig.savefig(W+"/results/fig2_artifact.png")
fig.savefig(W+"/results/fig2_artifact.pdf")
print("wrote fig2_artifact.png and .pdf")
for c in CLS:
    print("%-4s n=%-9d complete %.1f%%  mean all %.0f  mean complete %.0f  antifam %.3f%%"%(
        c,n[c],100.0*comp[c]["00"]/n[c],np.mean(L[c]),np.mean(Lc[c]),100.0*sp[c]/n[c]))
