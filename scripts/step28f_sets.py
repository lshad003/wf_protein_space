# Step 28f: calibration query sets for the ProstT5/Foldseek pilot (step28e).
# a. results/prostt5_calib_known_lenmatched.faa: 740 K-class (Pfam domain) supported reps from
#    results/fourway_classes_n88.tsv, completeness 00 and >= 35 aa (same pool as step28d), drawn with
#    random.Random(28) per length bin to match the bin counts of results/prostt5_pilot_eu740.faa
#    (<100, 100-199, 200-399, >=400 aa).
# b. results/prostt5_decoy_eu740rev.faa: the 740 pilot sequences reversed, headers prefixed "rev_".
# aa lengths from catalog/db/LsPS_AA_all.fasta.gz.fai; sequences streamed from the .gz. Refuses to overwrite.
import os, gzip, random, statistics
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OA=W+"/results/prostt5_calib_known_lenmatched.faa"; OB=W+"/results/prostt5_decoy_eu740rev.faa"
for o in (OA,OB):
    if os.path.exists(o): raise SystemExit("REFUSING: "+o+" exists")
BINS=(("<100",35,100),("100-199",100,200),("200-399",200,400),(">=400",400,10**9))
def binof(n):
    for b,lo,hi in BINS:
        if lo<=n<hi: return b
pil=[]; h=None
for l in open(W+"/results/prostt5_pilot_eu740.faa"):
    if l[0]==">": h=l[1:].strip()
    else: pil.append((h,l.strip()))
if len(pil)!=740: raise SystemExit("ABORT: pilot is not 740 sequences")
need={b:sum(1 for _,s in pil if binof(len(s))==b) for b,_,_ in BINS}
print("pilot bin counts:", need, flush=True)
aalen={}
for l in open(W+"/catalog/db/LsPS_AA_all.fasta.gz.fai"):
    f=l.split("\t",2); aalen[f[0]]=int(f[1])
pool={b:[] for b,_,_ in BINS}
with open(W+"/results/fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        f=l.split("\t")
        if f[1]=="K" and f[3]=="00" and aalen[f[0]]>=35: pool[binof(aalen[f[0]])].append(f[0])
print("K complete >= 35 aa pool by bin:", {b:len(v) for b,v in pool.items()}, flush=True)
rng=random.Random(28); pick=[]
for b,_,_ in BINS: pick+=rng.sample(sorted(pool[b]),need[b])
want=set(pick)
seq={}; h=None; s=[]
with gzip.open(W+"/catalog/db/LsPS_AA_all.fasta.gz","rt") as fh:
    for l in fh:
        if l[0]==">":
            if h in want: seq[h]="".join(s)
            h=l[1:].split()[0]; s=[]
        elif h in want: s.append(l.strip())
    if h in want: seq[h]="".join(s)
if set(seq)!=want: raise SystemExit("ABORT: %d sequences not found"%len(want-set(seq)))
if any(len(seq[x])!=aalen[x] for x in pick): raise SystemExit("ABORT: length disagrees with fai")
with open(OA,"w") as fh:
    for x in pick: fh.write(">%s\n%s\n"%(x,seq[x]))
with open(OB,"w") as fh:
    for x,sq in pil: fh.write(">rev_%s\n%s\n"%(x,sq[::-1]))
def dist(lab,L):
    L=sorted(L); q=lambda p:L[min(len(L)-1,int(p*len(L)))]
    print("%s\tn=%d\tmin=%d\tp25=%d\tmedian=%g\tp75=%d\tmax=%d\tmean=%.1f\tbins=%s"%(lab,len(L),L[0],q(.25),
        statistics.median(L),q(.75),L[-1],statistics.mean(L),[sum(1 for x in L if binof(x)==b) for b,_,_ in BINS]))
dist("EU740",[len(s) for _,s in pil]); dist("K740_lenmatched",[len(seq[x]) for x in pick])
dist("EU740rev",[len(s) for _,s in pil])
print("wrote", OA, len(pick)); print("wrote", OB, len(pil))
