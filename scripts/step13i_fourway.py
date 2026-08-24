# Step 13i: four-way split of supported representatives.
# K   = Pfam domain (--cut_ga)
# KWP = eggNOG assignment, no Pfam domain
# GU  = no annotation, but hits NCBI ClusteredNR 20260128
#       PRIMARY CRITERION: E <= 1e-10, qcovhsp >= 50, scovhsp >= 50
# EU  = no annotation and no nr hit at that criterion
import collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
E,QC,SC = 1e-10, 50.0, 50.0

hit=set()
for ln in open(W+"/results/step13_nrclust20260128_hits.tsv"):
    p=ln.rstrip("\n").split("\t")
    if float(p[5])<=E and float(p[7])>=QC and float(p[8])>=SC:
        hit.add(p[0])
print("nr hits at primary criterion:",len(hit),flush=True)

spur=set()
for ln in open(W+"/results/antifam_hits.tblout"):
    if ln[0]!="#": spur.add(ln.split()[0])

L={}
for ln in open(W+"/results/supported_cds_lengths.tsv"):
    g,l=ln.rstrip("\n").split("\t"); L[g]=int(l)
comp={}
for ln in open(W+"/results/completeness_reps.tsv"):
    g,c=ln.rstrip("\n").split("\t"); comp[g]=c
yr={}; prev={}; cnt={}
import gzip
with gzip.open(W+"/results/prevalence_primary.tsv.gz","rt") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        yr[p[0]]=(int(p[3])>0)+(int(p[4])>0)+(int(p[5])>0)
        prev[p[0]]=int(p[2]); cnt[p[0]]=int(p[6])

out=open(W+"/results/fourway_classes.tsv","w")
out.write("rep\tclass4\tlength\tcompleteness\tn_samples\tn_years\ttotal_counts\tantifam\n")
st=collections.defaultdict(lambda:{"n":0,"len":0,"cmp":0,"prev":0,"cnt":0,"sp":0,"y3":0})
with open(W+"/results/rep_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if p[1]!="1": continue
        g,k=p[0],p[5]
        c4 = k if k in ("K","KWP") else ("GU" if g in hit else "EU")
        s=st[c4]; s["n"]+=1; s["len"]+=L.get(g,0); s["prev"]+=prev.get(g,0)
        s["cnt"]+=cnt.get(g,0); s["sp"]+= (1 if g in spur else 0)
        if comp.get(g)=="00": s["cmp"]+=1
        if yr.get(g,0)==3: s["y3"]+=1
        out.write("%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n"%(
            g,c4,L.get(g,0),comp.get(g,"NA"),prev.get(g,0),yr.get(g,0),
            cnt.get(g,0),1 if g in spur else 0))
out.close()

T=2069453
tot=sum(st[c]["cnt"] for c in st)
print("\ncriterion for GU: E<=1e-10, qcov>=50, scov>=50, NCBI ClusteredNR 20260128")
print("\nclass   n          %cat    mean_bp  complete%  mean_prev  3yr%    %reads  antifam%")
for c in ["K","KWP","GU","EU"]:
    s=st[c]
    print("%-6s  %-10d %-7.2f %-8.0f %-10.1f %-10.1f %-7.1f %-7.2f %.3f"%(
        c, s["n"], 100.0*s["n"]/T, s["len"]/s["n"], 100.0*s["cmp"]/s["n"],
        s["prev"]/s["n"], 100.0*s["y3"]/s["n"], 100.0*s["cnt"]/tot,
        100.0*s["sp"]/s["n"]))
print("\ntotal:", sum(st[c]["n"] for c in st), "expect", T)
assert sum(st[c]["n"] for c in st)==T
print("wrote results/fourway_classes.tsv")
