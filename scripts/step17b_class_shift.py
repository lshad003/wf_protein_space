# Step 17b: does the U (unannotated) share of gene abundance differ by
# treatment, within cohort? Non-parametric, no distribution assumption.
# WF22: Kruskal-Wallis across 3 groups, then pairwise Mann-Whitney.
# WF23: Mann-Whitney, 5 vs 4.
# WF24: Kruskal-Wallis across 7 codes.
from scipy.stats import kruskal, mannwhitneyu
import collections, itertools
W="/bigdata/stajichlab/lshad003/wf_protein_space"
tr={}
with open(W+"/metadata/samples_89_treatment.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t"); tr[p[0]]=(p[1],p[2])
dat=collections.defaultdict(lambda: collections.defaultdict(list))
with open(W+"/results/dark_abundance_by_sample.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        c,t=tr[p[0]]
        dat[c][t].append((float(p[4]),float(p[5])))  # raw U%, rpk U%
for c in ["WF22","WF23","WF24"]:
    print("\n=== %s ==="%c)
    groups=sorted(dat[c])
    for lab,i in [("raw U%",0),("length-normalized U%",1)]:
        vals=[[x[i] for x in dat[c][g]] for g in groups]
        for g,v in zip(groups,vals):
            v2=sorted(v); med=v2[len(v2)//2]
            print("  %-14s %-12s n=%d median=%.2f range %.2f-%.2f"%(lab,g,len(v),med,v2[0],v2[-1]))
        if len(groups)==2:
            s,p=mannwhitneyu(vals[0],vals[1],alternative="two-sided")
            print("  --> Mann-Whitney U=%.1f p=%.4f"%(s,p))
        else:
            s,p=kruskal(*vals)
            print("  --> Kruskal-Wallis H=%.3f p=%.4f"%(s,p))
            if p<0.05:
                for a,b in itertools.combinations(range(len(groups)),2):
                    s2,p2=mannwhitneyu(vals[a],vals[b],alternative="two-sided")
                    print("      %s vs %s p=%.4f"%(groups[a],groups[b],p2))
        print("")
