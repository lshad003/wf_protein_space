#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 4 --mem 128gb --time=8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step9.%j.log
# Extracts Prodigal partial= flags for every protein, then reports ORF
# completeness for the catalog as a whole, for the 95% representatives, and
# split by cluster size. Headers were stripped of the flag during renaming,
# so the flags are read from the original per-sample aa.fa.gz files.
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import gzip, os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"

stems = (open(W+"/metadata/wf22_stems_44.txt").read().split()
       + open(W+"/metadata/new45_stems.txt").read().split())
print("stems:", len(stems), flush=True)

print("loading representative IDs ...", flush=True)
reps=set()
for line in open(W+"/catalog/db/LsPS_AA_95_rep.fasta"):
    if line[0]==">": reps.add(line[1:].strip())
print("reps:", len(reps), flush=True)

allc=collections.Counter()      # partial code, all proteins
repc={}                         # rep id -> partial code
for i,s in enumerate(stems,1):
    p = f"{NEW}/{s}.aa.fa.gz" if os.path.exists(f"{NEW}/{s}.aa.fa.gz") else f"{OLD}/{s}.aa.fa.gz"
    with gzip.open(p,"rt") as fh:
        for line in fh:
            if line[0]!=">": continue
            tok=line[1:].split(None,1)
            pid=s+"__"+tok[0]
            code="NA"
            if len(tok)>1:
                j=tok[1].find("partial=")
                if j>=0: code=tok[1][j+8:j+10]
            allc[code]+=1
            if pid in reps: repc[pid]=code
    if i%10==0: print(" ",i,"stems done", flush=True)

def show(title, counter, total):
    print(f"\n=== {title} (n={total:,}) ===")
    for k in sorted(counter):
        lab={"00":"complete","10":"no start","01":"no stop","11":"both missing"}.get(k,k)
        print(f"  {k} {lab:14} {counter[k]:>12,}  {100.0*counter[k]/total:6.2f}%")
    print(f"  COMPLETE FRACTION: {100.0*counter.get('00',0)/total:.2f}%")

show("all proteins", allc, sum(allc.values()))
rc=collections.Counter(repc.values())
show("95% representatives", rc, sum(rc.values()))

print("\nloading cluster sizes ...", flush=True)
size=collections.Counter()
for line in open(W+"/results/clusters_95.tsv"):
    size[line.split("\t",1)[0]]+=1

bins=[(1,1,"singleton"),(2,2,"2"),(3,5,"3-5"),(6,20,"6-20"),(21,10**9,"21+")]
print("\n=== completeness of representatives by cluster size ===")
print(f"{'bin':12} {'clusters':>10} {'complete':>10} {'pct':>7}")
out=[]
for lo,hi,lab in bins:
    n=c=0
    for r,code in repc.items():
        s=size.get(r,0)
        if lo<=s<=hi:
            n+=1
            if code=="00": c+=1
    print(f"{lab:12} {n:>10,} {c:>10,} {100.0*c/n if n else 0:6.2f}%")
    out.append((lab,n,c))

with open(W+"/results/orf_completeness.tsv","w") as fh:
    fh.write("category\tn\tcomplete\tpct_complete\n")
    fh.write(f"all_proteins\t{sum(allc.values())}\t{allc.get('00',0)}\t{100.0*allc.get('00',0)/sum(allc.values()):.2f}\n")
    fh.write(f"representatives\t{sum(rc.values())}\t{rc.get('00',0)}\t{100.0*rc.get('00',0)/sum(rc.values()):.2f}\n")
    for lab,n,c in out:
        fh.write(f"cluster_size_{lab}\t{n}\t{c}\t{100.0*c/n if n else 0:.2f}\n")
print("\nwrote", W+"/results/orf_completeness.tsv")
PYEOF
