# Step 26f2: recurrence of EU-only direct 30% families (step26f) across metagenomes and cohorts.
# Families: results/unk_direct30_eu_families.tsv. Members: results/unk_direct30_clusters.tsv
# (family rep, member). Stem = member ID before "__"; cohort from metadata/samples_88_treatment.tsv.
# 95% representative of each member from results/clusters_95_n88.tsv (rep, member), restricted to
# the 384,260 reps in results/unk_filtered_aplakidou.faa.
# Per family: protein records, distinct 95% reps, distinct metagenomes, distinct cohorts,
# distinct animals (stem prefix before "."; WF22 15 animals from 44 stems, WF23/WF24 one stem each).
# Output: results/unk_direct30_eu_recurrence.tsv; refuses to overwrite.
import os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/unk_direct30_eu_recurrence.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")
coh={}
with open(W+"/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh:
        f=l.rstrip("\n").split("\t"); coh[f[0]]=f[1]
print("stems:", len(coh), collections.Counter(coh.values()), flush=True)
eu={}
with open(W+"/results/unk_direct30_eu_families.tsv") as fh:
    next(fh)
    for l in fh:
        f,n,r=l.rstrip("\n").split("\t"); eu[f]=(int(n),int(r))
print("EU-only families:", len(eu), flush=True)
reps={l[1:].split()[0] for l in open(W+"/results/unk_filtered_aplakidou.faa") if l[0]==">"}
m2r={}
with open(W+"/results/clusters_95_n88.tsv") as fh:
    for l in fh:
        r,m=l.rstrip("\n").split("\t")
        if r in reps: m2r[m]=r
fam=collections.defaultdict(lambda: [0,set(),set(),set(),set()])
nostem=0
for l in open(W+"/results/unk_direct30_clusters.tsv"):
    f,m=l.rstrip("\n").split("\t")
    if f not in eu: continue
    s=m.split("__",1)[0]
    if s not in coh: nostem+=1; continue
    e=fam[f]; e[0]+=1; e[1].add(m2r[m]); e[2].add(s); e[3].add(coh[s]); e[4].add(s.split(".",1)[0])
if nostem: raise SystemExit("ABORT: %d members with a stem not in metadata"%nostem)
if len(fam)!=len(eu): raise SystemExit("ABORT: EU-only families missing members")
bad=sum(1 for f,(n,r) in eu.items() if (fam[f][0],len(fam[f][1]))!=(n,r))
if bad: raise SystemExit("ABORT: %d families disagree with unk_direct30_eu_families.tsv"%bad)
print("cross-check vs unk_direct30_eu_families.tsv (n_proteins, n_reps95): OK", flush=True)
rows=[(f,v[0],len(v[1]),len(v[2]),len(v[3]),len(v[4])) for f,v in fam.items()]
rows.sort(key=lambda x:(-x[3],-x[1],x[0]))
with open(OUT,"w") as fh:
    fh.write("family_rep\tn_proteins\tn_reps95\tn_metagenomes\tn_cohorts\tn_animals\n")
    for r in rows: fh.write("%s\t%d\t%d\t%d\t%d\t%d\n"%r)
crit=[(">= 3 protein records from >= 2 metagenomes", lambda p,r,m,c,a: p>=3 and m>=2),
      (">= 3 distinct 95% reps from >= 2 metagenomes", lambda p,r,m,c,a: r>=3 and m>=2),
      (">= 5 metagenomes", lambda p,r,m,c,a: m>=5),
      (">= 10 metagenomes", lambda p,r,m,c,a: m>=10),
      (">= 2 cohorts", lambda p,r,m,c,a: c>=2),
      ("all 3 cohorts", lambda p,r,m,c,a: c>=3),
      (">= 10 metagenomes AND >= 2 cohorts", lambda p,r,m,c,a: m>=10 and c>=2),
      (">= 10 animals", lambda p,r,m,c,a: a>=10),
      (">= 10 animals AND >= 2 cohorts", lambda p,r,m,c,a: a>=10 and c>=2)]
N=len(rows)
print("\ncriterion\tfamilies\tof\tpct")
for lab,fn in crit:
    k=sum(1 for _,p,r,m,c,a in rows if fn(p,r,m,c,a))
    print("%s\t%d\t%d\t%.2f"%(lab,k,N,100.0*k/N))
print("\ncross-tab (rows: metagenomes; cols: cohorts 1/2/3)")
for lab,lo,hi in (("< 10 metagenomes",0,10),(">= 10 metagenomes",10,10**9)):
    print(lab+"\t"+"\t".join(str(sum(1 for _,p,r,m,c,a in rows if lo<=m<hi and c==k)) for k in (1,2,3)))
print("\nNOTE: every input 95% rep is already a supported cluster (>= 3 members from >= 2")
print("metagenomes), so the protein-record criterion was expected to pass nearly every family.")
print("It does not: support was counted on all 95% members, but step26e1 clustered only complete")
print(">= 35 aa members, so a supported 95% cluster can fall below 3 records or 2 metagenomes here.")
print("wrote", OUT, "rows:", N)
