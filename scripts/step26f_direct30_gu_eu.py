# Step 26f: GU/EU composition of the direct 30% families (step26e1).
# Each member protein in results/unk_direct30_clusters.tsv takes the class4 label of
# its 95% representative (results/clusters_95_n88.tsv: rep, member; restricted to the
# 384,260 reps in results/unk_filtered_aplakidou.faa) from results/fourway_classes_n88.tsv.
# Family type: EU_only (all members EU), GU_only (all GU), mixed (both).
# Size criteria are stated two ways: protein records and distinct 95% reps.
# Output: results/unk_direct30_eu_families.tsv (one row per EU-only family); refuses to overwrite.
import os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/unk_direct30_eu_families.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")
reps={l[1:].split()[0] for l in open(W+"/results/unk_filtered_aplakidou.faa") if l[0]==">"}
cls={}
with open(W+"/results/fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        f=l.split("\t",2)
        if f[0] in reps: cls[f[0]]=f[1]
print("filtered reps:", len(reps), "| labelled:", len(cls), collections.Counter(cls.values()), flush=True)
if len(cls)!=len(reps): raise SystemExit("ABORT: filtered reps without a class label")
m2r={}
with open(W+"/results/clusters_95_n88.tsv") as fh:
    for l in fh:
        r,m=l.rstrip("\n").split("\t")
        if r in cls: m2r[m]=r
print("members of filtered reps:", len(m2r), flush=True)
fam=collections.defaultdict(lambda: [0,0,set()])  # n_eu, n_gu, 95% reps
other=0; missing=0; n=0
for l in open(W+"/results/unk_direct30_clusters.tsv"):
    f,m=l.rstrip("\n").split("\t"); n+=1
    r=m2r.get(m)
    if r is None: missing+=1; continue
    c=cls[r]; e=fam[f]
    if c=="EU": e[0]+=1
    elif c=="GU": e[1]+=1
    else: other+=1
    e[2].add(r)
print("proteins:", n, "| unmapped:", missing, "| non-GU/EU labels:", other, "| families:", len(fam), flush=True)
if missing or other: raise SystemExit("ABORT: unmapped members or non-GU/EU labels")
typ=collections.Counter(); size=collections.Counter()
rows=[]
for f,(ne,ng,rs) in fam.items():
    t="EU_only" if ng==0 else ("GU_only" if ne==0 else "mixed")
    typ[t]+=1
    if t=="EU_only":
        rows.append((f,ne,len(rs)))
        for k in (3,10):
            if ne>=k: size[("proteins",k)]+=1
            if len(rs)>=k: size[("reps95",k)]+=1
print("family types (all members labelled):", dict(typ), "| total:", sum(typ.values()))
print("EU_only families >=3 protein records:", size[("proteins",3)], "| >=10:", size[("proteins",10)])
print("EU_only families >=3 distinct 95% reps:", size[("reps95",3)], "| >=10:", size[("reps95",10)])
rows.sort(key=lambda x:(-x[1],x[0]))
with open(OUT,"w") as fh:
    fh.write("family_rep\tn_proteins\tn_reps95\n")
    for r in rows: fh.write("%s\t%d\t%d\n"%r)
print("wrote", OUT, "rows:", len(rows))
