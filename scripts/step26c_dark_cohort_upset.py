# Step 26c: cohort (year) sharing of the 8,464 unknown 50% families with >=3 members
# (results/unk_clusters_50_n88.tsv, members are 95% representatives). Read-only on _n88.
# Presence: a family is present in a cohort if any protein of its constituent 95%
# clusters (results/clusters_95_n88.tsv, 88 metagenomes) comes from that cohort.
# Sets reported for all 8,464 families and for families with >=10 members
# (members = 95% representatives, the same unit as the >=3 cut).
# Method rules: raw sharing is effort-confounded (WF22 44, WF23 9, WF24 35 samples),
# so sharing is also given at matched effort, 9 samples per cohort, 1,000 seeded
# draws. Run check: same cohort on different runs (WF24 M005990 vs M006342) against
# different cohorts on one run (WF23 vs WF24, both UCB_20250426_M005990), 9 vs 9.
# Outputs: results/dark_family_cohorts.tsv (per family),
# results/dark_family_cohorts_summary.tsv, figures/dark_family_cohorts_upset{,_ge10}.{pdf,png}
import collections, os, random, sys
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import wfstyle
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/dark_family_cohorts.tsv"
SUM=W+"/results/dark_family_cohorts_summary.tsv"
FIG=W+"/figures/dark_family_cohorts_upset"
for p in (OUT,SUM,FIG+".pdf",FIG+".png",FIG+"_ge10.pdf",FIG+"_ge10.png"):
    if os.path.exists(p): raise SystemExit("REFUSING: "+p+" exists")
COH=["WF22","WF23","WF24"]

meta={}; run={}
for line in open(W+"/results/batch_table_88.tsv"):
    f=line.rstrip("\n").split("\t")
    if f[0]=="stem": continue
    meta[f[0]]=f[1]; run[f[0]]=f[2]
stems=sorted(meta); idx={s:i for i,s in enumerate(stems)}
print("stems:", len(stems), collections.Counter(meta.values()), flush=True)

fam=collections.defaultdict(list)
for line in open(W+"/results/unk_clusters_50_n88.tsv"):
    f,r=line.rstrip("\n").split("\t"); fam[f].append(r)
fams=[f for f in fam if len(fam[f])>=3]
print("50% families total:", len(fam), "| >=3 members:", len(fams), flush=True)
if len(fams)!=8464: raise SystemExit("ABORT: expected 8,464 families")
rep2fam={r:f for f in fams for r in fam[f]}

mask=collections.Counter(); nprot=collections.Counter()
for line in open(W+"/results/clusters_95_n88.tsv"):
    r,m=line.rstrip("\n").split("\t")
    f=rep2fam.get(r)
    if f is None: continue
    mask[f]|=1<<idx[m.split("__",1)[0]]; nprot[f]+=1
nomem=[f for f in fams if mask[f]==0]
print("member proteins:", sum(nprot.values()), "| families with no n88 member:", len(nomem), flush=True)

cmask={c:sum(1<<idx[s] for s in stems if meta[s]==c) for c in COH}
def pattern(m, cm):
    p=[c for c in COH if m & cm[c]]
    return "+".join(p) if p else "none"
PATS=["WF22","WF23","WF24","WF22+WF23","WF22+WF24","WF23+WF24","WF22+WF23+WF24"]
ge10=set(f for f in fams if len(fam[f])>=10)

with open(OUT,"w") as fh:
    fh.write("family_rep\tn_reps_95\tn_proteins\tn_samples_WF22\tn_samples_WF23\tn_samples_WF24\tpattern\n")
    for f in sorted(fams, key=lambda x:(-len(fam[x]),x)):
        ns=[bin(mask[f]&cmask[c]).count("1") for c in COH]
        fh.write(f"{f}\t{len(fam[f])}\t{nprot[f]}\t"+"\t".join(map(str,ns))+f"\t{pattern(mask[f],cmask)}\n")
print("wrote", OUT, flush=True)

rows=[]
def raw(sub, lab):
    c=collections.Counter(pattern(mask[f],cmask) for f in sub)
    print(f"\nRAW, all samples ({lab}, n={len(sub):,})")
    for p in PATS+["none"]:
        if p=="none" and not c[p]: continue
        print(f"  {p:16} {c[p]:>7,}  {100.0*c[p]/len(sub):6.2f}%")
        rows.append((lab,"raw_all_samples",p,c[p],f"{100.0*c[p]/len(sub):.2f}","",""))
    return c
craw={}
for lab,sub in (("all_8464",fams),("ge10_members",[f for f in fams if f in ge10])):
    craw[lab]=raw(sub,lab)

def draws(groups, sub, n=9, ndraw=1000, seed=26):
    rng=random.Random(seed)
    pool={g:[s for s in stems if groups[s]==g] for g in set(groups.values()) if g}
    res=collections.defaultdict(list)
    names=sorted(pool)
    for _ in range(ndraw):
        cm={g:sum(1<<idx[s] for s in rng.sample(pool[g],n)) for g in names}
        c=collections.Counter()
        for f in sub:
            p=[g for g in names if mask[f]&cm[g]]
            c["+".join(p) if p else "none"]+=1
        for k in set(list(c)+["none"]): res[k].append(c[k])
    return res, names
def show(res, sub, lab, method, order):
    print(f"\n{method} ({lab}, n={len(sub):,}), mean [2.5%, 97.5%] over 1,000 draws")
    for p in order:
        v=sorted(res.get(p,[0]*1000)); mu=sum(v)/len(v)
        lo,hi=v[int(0.025*len(v))],v[int(0.975*len(v))-1]
        print(f"  {p:28} {mu:>9.1f} [{lo:,}, {hi:,}]  {100.0*mu/len(sub):6.2f}%")
        rows.append((lab,method,p,f"{mu:.1f}",f"{100.0*mu/len(sub):.2f}",lo,hi))

for lab,sub in (("all_8464",fams),("ge10_members",[f for f in fams if f in ge10])):
    r,_=draws(meta,sub)
    show(r,sub,lab,"matched_9_per_cohort",PATS+["none"])
    g={s:("WF24_M005990" if run[s]=="UCB_20250426_M005990" else "WF24_M006342") if meta[s]=="WF24" else None for s in stems}
    r,n=draws(g,sub)
    show(r,sub,lab,"run_check_same_cohort_diff_run",[n[0],n[1],"+".join(n),"none"])
    g={s:(meta[s]+"_M005990" if run[s]=="UCB_20250426_M005990" else None) for s in stems}
    r,n=draws(g,sub)
    show(r,sub,lab,"run_check_diff_cohort_same_run",[n[0],n[1],"+".join(n),"none"])

with open(SUM,"w") as fh:
    fh.write("family_set\tmethod\tpattern\tn_families\tpct\tci_low\tci_high\n")
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
print("\nwrote", SUM)

def upset(sub, path, title):
    pats=[pattern(mask[f],cmask).split("+") for f in sub]
    data=from_memberships(pats, data=[1]*len(pats)).groupby(level=list(range(3))).count()
    fig=plt.figure(figsize=wfstyle.size(120,80))
    u=UpSet(data, sort_by="cardinality", show_counts="{:,}", facecolor="#333333", element_size=None)
    for c in COH: u.style_categories(c, bar_facecolor=wfstyle.COH[c], bar_edgecolor="#333333")
    u.plot(fig=fig)
    fig.suptitle(title, fontsize=9)
    for ext in ("pdf","png"): fig.savefig(path+"."+ext)
    plt.close(fig); print("wrote", path+".pdf/.png")
upset(fams, FIG, f"Unknown 50% families, >=3 members (n={len(fams):,}), all 88 samples")
upset([f for f in fams if f in ge10], FIG+"_ge10", f"Unknown 50% families, >=10 members (n={len(ge10):,}), all 88 samples")
print("DONE step26c")
