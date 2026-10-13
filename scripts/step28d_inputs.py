# Step 28d: input FASTAs for the ProstT5/Foldseek pilot.
# 1. results/prostt5_pilot_eu740.faa: the linclust family representative of each EU-only 30%
#    family with >= 3 distinct 95% reps from >= 2 metagenomes (results/unk_direct30_eu_recurrence.tsv).
# 2. results/prostt5_calib_known1000.faa: 1,000 random (seed 28) K-class (Pfam domain) supported reps
#    from results/fourway_classes_n88.tsv with completeness 00 (complete ORF) and >= 35 aa.
# Lengths in aa from catalog/db/LsPS_AA_all.fasta.gz.fai (the length column of the four-way table is
# CDS bp). Sequences streamed from catalog/db/LsPS_AA_all.fasta.gz. Refuses to overwrite.
import os, gzip, random, statistics
W="/bigdata/stajichlab/lshad003/wf_protein_space"
O1=W+"/results/prostt5_pilot_eu740.faa"; O2=W+"/results/prostt5_calib_known1000.faa"
for o in (O1,O2):
    if os.path.exists(o): raise SystemExit("REFUSING: "+o+" exists")
aalen={}
for l in open(W+"/catalog/db/LsPS_AA_all.fasta.gz.fai"):
    f=l.split("\t",2); aalen[f[0]]=int(f[1])
print("proteins in fai:", len(aalen), flush=True)
eu=[]
with open(W+"/results/unk_direct30_eu_recurrence.tsv") as fh:
    next(fh)
    for l in fh:
        f=l.rstrip("\n").split("\t")
        if int(f[2])>=3 and int(f[3])>=2: eu.append(f[0])
print("EU families >= 3 distinct 95% reps from >= 2 metagenomes:", len(eu), flush=True)
if len(eu)!=740: raise SystemExit("ABORT: expected 740")
short=[x for x in eu if aalen[x]<35]
if short: raise SystemExit("ABORT: %d EU family reps < 35 aa"%len(short))
pool=[]; nk=0
with open(W+"/results/fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        f=l.split("\t")
        if f[1]!="K": continue
        nk+=1
        if f[3]=="00" and aalen[f[0]]>=35: pool.append(f[0])
print("K reps:", nk, "| complete and >= 35 aa:", len(pool), flush=True)
cal=random.Random(28).sample(sorted(pool),1000)
want=set(eu)|set(cal)
if len(want)!=1740: raise SystemExit("ABORT: overlap between sets")
seq={}; h=None; s=[]
with gzip.open(W+"/catalog/db/LsPS_AA_all.fasta.gz","rt") as fh:
    for l in fh:
        if l[0]==">":
            if h in want: seq[h]="".join(s)
            h=l[1:].split()[0]; s=[]
        elif h in want: s.append(l.strip())
    if h in want: seq[h]="".join(s)
miss=want-set(seq)
if miss: raise SystemExit("ABORT: %d sequences not found"%len(miss))
bad=[x for x in want if len(seq[x])!=aalen[x]]
if bad: raise SystemExit("ABORT: %d sequence lengths disagree with fai"%len(bad))
def dist(lab,ids):
    L=sorted(len(seq[x]) for x in ids)
    q=lambda p:L[min(len(L)-1,int(p*len(L)))]
    print("%s\tn=%d\tmin=%d\tp10=%d\tp25=%d\tmedian=%g\tp75=%d\tp90=%d\tmax=%d\tmean=%.1f"%(
        lab,len(L),L[0],q(.1),q(.25),statistics.median(L),q(.75),q(.9),L[-1],statistics.mean(L)))
    print("%s\t<100: %d\t100-199: %d\t200-399: %d\t>=400: %d"%(lab,
        sum(x<100 for x in L),sum(100<=x<200 for x in L),sum(200<=x<400 for x in L),sum(x>=400 for x in L)))
for o,ids in ((O1,eu),(O2,cal)):
    with open(o,"w") as fh:
        for x in ids: fh.write(">%s\n%s\n"%(x,seq[x]))
    print("wrote", o, len(ids))
dist("EU740",eu); dist("K1000",cal)
