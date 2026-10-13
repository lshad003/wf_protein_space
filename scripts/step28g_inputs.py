# Step 28g: main ProstT5/Foldseek query set.
# results/prostt5_main_eu15580.faa: the linclust family representative (as step28d) of each EU-only
# 30% family with >= 10 animals AND >= 2 cohorts in results/unk_direct30_eu_recurrence.tsv.
# Aborts unless 15,580 families. Sequences streamed from catalog/db/LsPS_AA_all.fasta.gz, lengths
# checked against its .fai. Refuses to overwrite.
import os, gzip, statistics
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/prostt5_main_eu15580.faa"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")
ids=[]
with open(W+"/results/unk_direct30_eu_recurrence.tsv") as fh:
    hdr=next(fh).rstrip("\n").split("\t")
    ic,ia=hdr.index("n_cohorts"),hdr.index("n_animals")
    for l in fh:
        f=l.rstrip("\n").split("\t")
        if int(f[ia])>=10 and int(f[ic])>=2: ids.append(f[0])
print("EU-only families >= 10 animals AND >= 2 cohorts:", len(ids), flush=True)
if len(ids)!=15580: raise SystemExit("ABORT: expected 15,580")
want=set(ids)
aalen={}
for l in open(W+"/catalog/db/LsPS_AA_all.fasta.gz.fai"):
    f=l.split("\t",2)
    if f[0] in want: aalen[f[0]]=int(f[1])
seq={}; h=None; s=[]
with gzip.open(W+"/catalog/db/LsPS_AA_all.fasta.gz","rt") as fh:
    for l in fh:
        if l[0]==">":
            if h in want: seq[h]="".join(s)
            h=l[1:].split()[0]; s=[]
        elif h in want: s.append(l.strip())
    if h in want: seq[h]="".join(s)
if set(seq)!=want: raise SystemExit("ABORT: %d sequences not found"%len(want-set(seq)))
if any(len(seq[x])!=aalen[x] for x in ids): raise SystemExit("ABORT: length disagrees with fai")
if any(len(seq[x])<35 for x in ids): raise SystemExit("ABORT: sequence < 35 aa")
with open(OUT,"w") as fh:
    for x in ids: fh.write(">%s\n%s\n"%(x,seq[x]))
L=sorted(len(seq[x]) for x in ids); q=lambda p:L[min(len(L)-1,int(p*len(L)))]
print("lengths aa: n=%d min=%d p25=%d median=%g p75=%d max=%d mean=%.1f"%(len(L),L[0],q(.25),statistics.median(L),q(.75),L[-1],statistics.mean(L)))
print("bins <100: %d 100-199: %d 200-399: %d >=400: %d"%(sum(x<100 for x in L),sum(100<=x<200 for x in L),sum(200<=x<400 for x in L),sum(x>=400 for x in L)))
print("wrote", OUT, len(ids))
