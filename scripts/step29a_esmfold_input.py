# Step 29a: ESMFold input. One family rep each from the 15,580 EU-only families (results/prostt5_main_eu15580.faa):
#  HIT:   best E over PDB and CATH50 <= 1e-3 (results/prostt5_main/eu15580_vs_{PDB,CATH50}.tsv); expect 218.
#  NOHIT: no hit at E <= 1e-1 in either DB; top 200 by n_animals, then n_metagenomes (desc), then rep ID,
#         from results/unk_direct30_eu_recurrence.tsv.
#  NOHIT_LONG: same no-hit pool, length >= 100 aa, not already in NOHIT; top 200, same ranking.
# Terminal '*' removed (ESMFold alphabet has no stop). Length = residues without '*'.
# Outputs: results/esmfold_input.faa, results/esmfold_input_labels.tsv; refuses to overwrite.
import os
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/esmfold_input.faa"; LAB=W+"/results/esmfold_input_labels.tsv"
for p in (OUT,LAB):
    if os.path.exists(p): raise SystemExit("REFUSING: "+p+" exists")
seq={}; h=None
for l in open(W+"/results/prostt5_main_eu15580.faa"):
    if l[0]==">": h=l[1:].split()[0]; seq[h]=""
    else: seq[h]+=l.strip()
if len(seq)!=15580: raise SystemExit("ABORT: query set is not 15,580")
best={}
for db in ("PDB","CATH50"):
    for l in open(W+"/results/prostt5_main/eu15580_vs_%s.tsv"%db):
        f=l.split("\t")
        if f[0] not in seq: raise SystemExit("ABORT: hit query not in query set")
        e=float(f[2]); best[f[0]]=min(e,best.get(f[0],e))
hit=sorted(x for x,e in best.items() if e<=1e-3)
if len(hit)!=218: raise SystemExit("ABORT: expected 218 HIT, got %d"%len(hit))
rec={}
with open(W+"/results/unk_direct30_eu_recurrence.tsv") as fh:
    hdr=next(fh).rstrip("\n").split("\t")
    ia,im=hdr.index("n_animals"),hdr.index("n_metagenomes")
    for l in fh:
        f=l.rstrip("\n").split("\t")
        if f[0] in seq: rec[f[0]]=(int(f[ia]),int(f[im]))
if set(rec)!=set(seq): raise SystemExit("ABORT: %d query families missing from recurrence table"%len(set(seq)-set(rec)))
pool=[x for x in seq if best.get(x,9e9)>1e-1]
pool.sort(key=lambda x:(-rec[x][0],-rec[x][1],x))
nohit=pool[:200]
print("families with no hit at E <= 1e-1 (either DB): %d"%len(pool))
L=lambda x:len(seq[x].rstrip("*"))
lpool=[x for x in pool if L(x)>=100 and x not in set(nohit)]
nohit_long=lpool[:200]
print("no-hit families >= 100 aa not in NOHIT: %d"%len(lpool))
print("NOHIT rank-200 cutoff: n_animals %d, n_metagenomes %d; rank 201: %s"%(rec[nohit[-1]]+(("n_animals %d, n_metagenomes %d"%rec[pool[200]]) if len(pool)>200 else "none",)))
print("NOHIT_LONG rank-200 cutoff: n_animals %d, n_metagenomes %d; rank 201: %s"%(rec[nohit_long[-1]]+(("n_animals %d, n_metagenomes %d"%rec[lpool[200]]) if len(lpool)>200 else "none",)))
rows=[(x,"HIT") for x in hit]+[(x,"NOHIT") for x in nohit]+[(x,"NOHIT_LONG") for x in nohit_long]
if len(set(r[0] for r in rows))!=len(rows): raise SystemExit("ABORT: duplicate rep across groups")
with open(OUT,"w") as fa, open(LAB,"w") as fl:
    fl.write("rep\tgroup\tlength\tn_animals\n")
    for x,g in rows:
        s=seq[x].rstrip("*")
        if "*" in s: raise SystemExit("ABORT: internal stop in "+x)
        fa.write(">%s\n%s\n"%(x,s)); fl.write("%s\t%s\t%d\t%d\n"%(x,g,len(s),rec[x][0]))
import statistics
for g,ids in (("HIT",hit),("NOHIT",nohit),("NOHIT_LONG",nohit_long)):
    L=[len(seq[x].rstrip("*")) for x in ids]; A=[rec[x][0] for x in ids]
    print("%s n=%d length median %g min %d max %d sum %d; n_animals median %g min %d max %d"%(g,len(ids),statistics.median(L),min(L),max(L),sum(L),statistics.median(A),min(A),max(A)))
print("wrote", OUT, len(rows)); print("wrote", LAB)
