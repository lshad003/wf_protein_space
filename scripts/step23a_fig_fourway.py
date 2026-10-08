# Figure: the four-way split and how the classes differ.
# Panels: A class sizes, B mean length, C completeness, D read share.
import collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
W="/bigdata/stajichlab/lshad003/wf_protein_space"
C={"K":"#0072B2","KWP":"#56B4E9","GU":"#E69F00","EU":"#D55E00"}
O=["K","KWP","GU","EU"]

st=collections.defaultdict(lambda:{"n":0,"len":0,"cmp":0,"cnt":0,"y3":0})
with open(W+"/results/fourway_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        s=st[p[1]]; s["n"]+=1; s["len"]+=int(p[2]); s["cnt"]+=int(p[6])
        if p[3]=="00": s["cmp"]+=1
        if p[5]=="3": s["y3"]+=1
tot=sum(st[c]["cnt"] for c in O); T=sum(st[c]["n"] for c in O)
print("total:",T)

fig,ax=plt.subplots(2,2,figsize=(9,7))
def bar(a,vals,ylab,title,fmt="%.1f"):
    a.bar(range(4),vals,color=[C[c] for c in O])
    a.set_xticks(range(4)); a.set_xticklabels(O)
    a.set_ylabel(ylab); a.set_title(title,fontsize=10)
    for i,v in enumerate(vals): a.text(i,v,fmt%v,ha="center",va="bottom",fontsize=8)
    for s in ("top","right"): a.spines[s].set_visible(False)

bar(ax[0,0],[100.0*st[c]["n"]/T for c in O],"Percent of catalog",
    "A  Class size (n = %s genes)"%format(T,","))
bar(ax[0,1],[st[c]["len"]/st[c]["n"] for c in O],"Mean CDS length (bp)",
    "B  Gene length","%.0f")
bar(ax[1,0],[100.0*st[c]["cmp"]/st[c]["n"] for c in O],"Percent complete ORFs",
    "C  Reading frame completeness")
bar(ax[1,1],[100.0*st[c]["cnt"]/tot for c in O],"Percent of mapped reads",
    "D  Share of gene abundance","%.2f")
fig.suptitle("Four-way classification of the support-filtered catalog",fontsize=11)
fig.tight_layout()
fig.savefig(W+"/results/fig_fourway.png",dpi=300)
print("wrote results/fig_fourway.png")
for c in O:
    s=st[c]
    print("%-4s n=%-9d %.2f%%  len %.0f  complete %.1f%%  reads %.2f%%  3yr %.1f%%"%(
        c,s["n"],100.0*s["n"]/T,s["len"]/s["n"],100.0*s["cmp"]/s["n"],
        100.0*s["cnt"]/tot,100.0*s["y3"]/s["n"]))
