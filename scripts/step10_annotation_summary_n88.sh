#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 4 --mem 200gb --time=8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step10_n88.%j.log
# n88 copy: representatives are the 95% clusters kept in results/reps_95_n88.txt,
# cluster sizes from results/clusters_95_n88.tsv, all-protein totals exclude
# UHM586.41010. Flags are still read from all 89 protein files, because a kept
# representative can be a UHM586.41010 sequence.
# eggNOG hits are restricted to that representative set; families from clusters_50_n88.tsv.
# Merges the 619 eggNOG chunk outputs, reports the annotated fraction of the
# 95% representatives, and breaks it down by cluster size, by ORF completeness,
# and rolled up to the 50% family tier. Completeness flags are re-read from the
# original per-sample protein files, since renaming removed them.
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import glob, gzip, os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"

files=sorted(glob.glob(W+"/results/eggnog2/*.emapper.annotations"))
print("annotation files:", len(files), flush=True)
annot=set()
empty=0
for f in files:
    n=0
    for line in open(f):
        if line.startswith("#"): continue
        annot.add(line.split("\t",1)[0]); n+=1
    if n==0: empty+=1
print("files with zero rows:", empty)
print("annotated representatives, 89 catalog:", len(annot), flush=True)

reps=[]
for line in open(W+"/results/reps_95_n88.txt"):
    reps.append(line.strip())
NR=len(reps)
annot&=set(reps)
print("annotated representatives, n88:", len(annot), flush=True)
print("total representatives:", NR)
print()
print("=== HEADLINE ===")
print(f"annotated:   {len(annot):>10,}  {100.0*len(annot)/NR:6.2f}%")
print(f"unannotated: {NR-len(annot):>10,}  {100.0*(NR-len(annot))/NR:6.2f}%")

print("\nloading cluster sizes ...", flush=True)
size=collections.Counter()
for line in open(W+"/results/clusters_95_n88.tsv"):
    size[line.split("\t",1)[0]]+=1

print("\n=== annotation rate by cluster size ===")
print(f"{'bin':12} {'clusters':>10} {'annotated':>10} {'pct':>7}")
rows=[]
for lo,hi,lab in [(1,1,"singleton"),(2,2,"2"),(3,5,"3-5"),(6,20,"6-20"),(21,10**9,"21+")]:
    n=a=0
    for r in reps:
        s=size.get(r,0)
        if lo<=s<=hi:
            n+=1
            if r in annot: a+=1
    print(f"{lab:12} {n:>10,} {a:>10,} {100.0*a/n if n else 0:6.2f}%")
    rows.append(("size_"+lab,n,a))

print("\nreading ORF completeness flags ...", flush=True)
OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
repset=set(reps)
comp={}
stems=(open(W+"/metadata/wf22_stems_44.txt").read().split()
     + open(W+"/metadata/new45_stems.txt").read().split())
for s in stems:
    p = f"{NEW}/{s}.aa.fa.gz" if os.path.exists(f"{NEW}/{s}.aa.fa.gz") else f"{OLD}/{s}.aa.fa.gz"
    with gzip.open(p,"rt") as fh:
        for line in fh:
            if line[0]!=">": continue
            tok=line[1:].split(None,1)
            pid=s+"__"+tok[0]
            if pid in repset and len(tok)>1:
                j=tok[1].find("partial=")
                if j>=0: comp[pid]=tok[1][j+8:j+10]

print("\n=== annotation rate by ORF completeness ===")
for code,lab in [("00","complete"),("10","no start"),("01","no stop"),("11","both missing")]:
    ids=[r for r in reps if comp.get(r)==code]
    a=sum(1 for r in ids if r in annot)
    if ids:
        print(f"  {lab:14} {len(ids):>10,} {a:>10,}  {100.0*a/len(ids):6.2f}%")
        rows.append(("orf_"+lab.replace(" ","_"),len(ids),a))

print("\nrolling up to the 50% family tier ...", flush=True)
fam=collections.defaultdict(list)
for line in open(W+"/results/clusters_50_n88.tsv"):
    rep,mem=line.rstrip("\n").split("\t")[:2]
    if mem in repset: fam[rep].append(mem)
nf=len(fam)
fa=sum(1 for r,v in fam.items() if any(m in annot for m in v))
print(f"\n=== 50% family tier ===")
print(f"families containing a 95% representative: {nf:,}")
print(f"  with at least one annotated member: {fa:,}  {100.0*fa/nf:6.2f}%")
print(f"  with no annotated member:           {nf-fa:,}  {100.0*(nf-fa)/nf:6.2f}%")
rows.append(("family_tier_50",nf,fa))

with open(W+"/results/annotation_summary_n88.tsv","w") as fh:
    fh.write("category\tn\tannotated\tpct_annotated\n")
    fh.write(f"all_representatives\t{NR}\t{len(annot)}\t{100.0*len(annot)/NR:.2f}\n")
    for lab,n,a in rows:
        fh.write(f"{lab}\t{n}\t{a}\t{100.0*a/n if n else 0:.2f}\n")
print("\nwrote", W+"/results/annotation_summary_n88.tsv")
PYEOF
