# Figure: maternal origin structures the catalog, fungal exposure does not.
# Numbers are transcribed from verified test output, source noted per panel.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
W="/bigdata/stajichlab/lshad003/wf_protein_space"

fig,ax=plt.subplots(2,2,figsize=(9.5,7))

# A: blocking changes the answer (logs/step17g.log)
lab=["STP1710.7","STP1717.1"]; nb=[4993,7275]; bl=[3,155]
x=np.arange(2); w=0.38
ax[0,0].bar(x-w/2,nb,w,color="#999999",label="samples treated as independent")
ax[0,0].bar(x+w/2,bl,w,color="#0072B2",label="blocked on animal")
ax[0,0].set_yscale("log"); ax[0,0].set_ylim(0.5,2e4)
ax[0,0].set_xticks(x); ax[0,0].set_xticklabels(lab)
ax[0,0].set_ylabel("Genes with padj < 0.05")
ax[0,0].legend(frameon=False,fontsize=7,loc="upper left")
ax[0,0].set_title("A  WF22: 44 samples from 15 animals",fontsize=10)
for i,(a,b) in enumerate(zip(nb,bl)):
    ax[0,0].text(i-w/2,a,str(a),ha="center",va="bottom",fontsize=8)
    ax[0,0].text(i+w/2,b,str(b),ha="center",va="bottom",fontsize=8)

# B: variance explained, animal centroids (logs/step22b.log)
n2=["Egg mass","Month","Treatment"]; v2=[51.3,12.7,6.9]
p2=["p = 0.001","p = 0.001","p = 0.657"]
c2=["#009E73","#56B4E9","#D55E00"]
ax[0,1].barh(range(3),v2,color=c2)
ax[0,1].set_yticks(range(3)); ax[0,1].set_yticklabels(n2)
ax[0,1].invert_yaxis(); ax[0,1].set_xlabel("Percent of variation explained (R2)")
ax[0,1].set_xlim(0,62)
ax[0,1].set_title("B  WF22 community composition",fontsize=10)
for i,(v,p) in enumerate(zip(v2,p2)):
    ax[0,1].text(v+1,i,"%.1f%%  %s"%(v,p),va="center",fontsize=8)

# C: every treatment contrast (logs/step17g, step17h)
nm=["WF22 STP1710.7","WF22 STP1717.1","WF22 EM3 only","WF23 UHM520",
    "WF24 code1","WF24 code2","WF24 code3","WF24 code5","WF24 code6","WF24 code7"]
vv=[3,155,0,0,0,0,1,0,0,0]
ax[1,0].barh(range(len(nm)),[max(v,0.25) for v in vv],color="#D55E00")
ax[1,0].set_yticks(range(len(nm))); ax[1,0].set_yticklabels(nm,fontsize=7.5)
ax[1,0].invert_yaxis(); ax[1,0].set_xlabel("Genes with padj < 0.05")
ax[1,0].set_xlim(0,180)
ax[1,0].set_title("C  All treatment contrasts, blocked where needed",fontsize=10)
for i,v in enumerate(vv): ax[1,0].text(max(v,0.25)+3,i,str(v),va="center",fontsize=8)

# D: same test, four datasets (logs/step22b, step22c, step22d)
nm4=["Abundance\n(WF22)","Presence\n(WF22)","Carbohydrate\nsubset (WF22)","Abundance\n(WF24)"]
egg=[51.3,36.3,67.3,np.nan]; trt=[6.9,10.6,5.3,14.9]
x=np.arange(4); w=0.38
ax[1,1].bar(x-w/2,egg,w,color="#009E73",label="egg mass")
ax[1,1].bar(x+w/2,trt,w,color="#D55E00",label="treatment")
ax[1,1].set_xticks(x); ax[1,1].set_xticklabels(nm4,fontsize=7.5)
ax[1,1].set_ylabel("Percent of variation explained (R2)")
ax[1,1].legend(frameon=False,fontsize=8)
ax[1,1].set_title("D  Treatment is null under every test",fontsize=10)
ax[1,1].text(3-w/2,2,"n/a",ha="center",fontsize=7,color="#666666")

for a in ax.flat:
    for s in ("top","right"): a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(W+"/results/fig_treatment_vs_eggmass.png",dpi=300)
print("wrote results/fig_treatment_vs_eggmass.png")
