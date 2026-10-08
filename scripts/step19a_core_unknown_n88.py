# Step 19a n88: copy of step19a_core_unknown.py on the 88-metagenome catalog:
# classes from results/rep_classes_n88.tsv, prevalence from
# results/prevalence_primary_n88.tsv.gz, AntiFam from results/antifam_hits_n88.tblout.
# Step 19a: characterize unknown genes by how many years they appear in.
# Compares the 3-year core against 2-year and 1-year unknowns on
# length, ORF completeness, prevalence, and abundance share.
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
comp={}
for ln in open(W+"/results/completeness_reps.tsv"):
    g,c=ln.rstrip("\n").split("\t"); comp[g]=c
spur=set()
for ln in open(W+"/results/antifam_hits_n88.tblout"):
    if ln[0]!="#": spur.add(ln.split()[0])
print("antifam flagged:",len(spur),flush=True)

# years per gene, plus prevalence and total counts, from prevalence table
out=gzip.open(W+"/results/core_unknown_by_years_n88.tsv.gz","wt")
out.write("gene\tclass\tn_years\tn_samples\ttotal_counts\tlength\tcompleteness\tantifam\n")
stat=collections.defaultdict(lambda: {"n":0,"len":0,"cmpl":0,"prev":0,"cnt":0,"spur":0,"clen":0,"cn":0})
with gzip.open(W+"/results/prevalence_primary_n88.tsv.gz","rt") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        g,k=p[0],p[1]
        ns=int(p[2]); y=(int(p[3])>0)+(int(p[4])>0)+(int(p[5])>0)
        tot=int(p[6]); ln_bp=L.get(g,0); c=comp.get(g,"NA")
        sp=1 if g in spur else 0
        out.write("%s\t%s\t%d\t%d\t%d\t%d\t%s\t%d\n"%(g,k,y,ns,tot,ln_bp,c,sp))
        if k!="U": continue
        s=stat[y]
        s["n"]+=1; s["len"]+=ln_bp; s["prev"]+=ns; s["cnt"]+=tot; s["spur"]+=sp
        if c=="00":
            s["cmpl"]+=1; s["clen"]+=ln_bp; s["cn"]+=1
out.close()

print("\n=== unknown (U) genes by number of years detected ===")
print("years    n        mean_bp  complete%  mean_prev  spurious%  complete_mean_bp")
for y in [3,2,1,0]:
    s=stat[y]
    if not s["n"]: continue
    print("%-8d %-8d %-8.0f %-10.1f %-10.1f %-10.3f %.0f"%(
        y, s["n"], s["len"]/s["n"], 100.0*s["cmpl"]/s["n"], s["prev"]/s["n"],
        100.0*s["spur"]/s["n"], s["clen"]/s["cn"] if s["cn"] else 0))
tc=sum(stat[y]["cnt"] for y in stat)
print("\n=== share of unknown-gene reads by year class ===")
for y in [3,2,1,0]:
    if stat[y]["n"]: print("  %d years: %.2f%% of U reads"%(y,100.0*stat[y]["cnt"]/tc))
print("\nwrote results/core_unknown_by_years_n88.tsv.gz")
