#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 4 --mem 220gb --time=8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step11.%j.log
# Two questions. First, whether the annotation rate of contig-spanning proteins
# (partial=11) is explained by cluster size rather than by completeness, tested
# as an interaction. Second, the unannotated fraction at the 50% family tier
# after support filtering, which is the figure comparable to published catalogs.
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import glob, gzip, os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"

annot=set()
for f in sorted(glob.glob(W+"/results/eggnog2/*.emapper.annotations")):
    for line in open(f):
        if not line.startswith("#"): annot.add(line.split("\t",1)[0])
print("annotated reps:", len(annot), flush=True)

repset=set()
for line in open(W+"/catalog/db/LsPS_AA_95_rep.fasta"):
    if line[0]==">": repset.add(line[1:].strip())

size=collections.Counter()
for line in open(W+"/results/clusters_95.tsv"):
    size[line.split("\t",1)[0]]+=1
print("cluster sizes loaded", flush=True)

OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
stems=(open(W+"/metadata/wf22_stems_44.txt").read().split()
     + open(W+"/metadata/new45_stems.txt").read().split())
comp={}
for s in stems:
    p=f"{NEW}/{s}.aa.fa.gz" if os.path.exists(f"{NEW}/{s}.aa.fa.gz") else f"{OLD}/{s}.aa.fa.gz"
    with gzip.open(p,"rt") as fh:
        for line in fh:
            if line[0]!=">": continue
            tok=line[1:].split(None,1)
            pid=s+"__"+tok[0]
            if pid in repset and len(tok)>1:
                j=tok[1].find("partial=")
                if j>=0: comp[pid]=tok[1][j+8:j+10]
print("completeness flags loaded", flush=True)

bins=[(1,1,"singleton"),(2,2,"2"),(3,5,"3-5"),(6,20,"6-20"),(21,10**9,"21+")]
print("\n=== cluster-size distribution of each completeness class ===")
print(f"{'flag':14}" + "".join(f"{lab:>12}" for _,_,lab in bins))
for code,lab in [("00","complete"),("10","no start"),("01","no stop"),("11","spans contig")]:
    ids=[r for r in repset if comp.get(r)==code]
    row=[]
    for lo,hi,_ in bins:
        row.append(sum(1 for r in ids if lo<=size.get(r,0)<=hi))
    tot=len(ids)
    print(f"{lab:14}" + "".join(f"{100.0*x/tot:11.1f}%" for x in row))

print("\n=== annotation rate, completeness WITHIN cluster size ===")
print(f"{'size bin':12}" + "".join(f"{l:>14}" for l in ["complete","no start","no stop","spans contig"]))
for lo,hi,lab in bins:
    cells=[]
    for code in ("00","10","01","11"):
        ids=[r for r in repset if comp.get(r)==code and lo<=size.get(r,0)<=hi]
        a=sum(1 for r in ids if r in annot)
        cells.append(f"{100.0*a/len(ids):12.1f}%" if ids else "          na")
    print(f"{lab:12}" + "".join(cells))

print("\nloading 50% families ...", flush=True)
fam_n=collections.Counter(); fam_s=collections.defaultdict(set); fam_rep=collections.defaultdict(list)
for line in open(W+"/results/clusters_50.tsv"):
    rep,mem=line.rstrip("\n").split("\t")[:2]
    fam_n[rep]+=1
    fam_s[rep].add(mem.split("__")[0])
    if mem in repset: fam_rep[rep].append(mem)

def report(title, keep):
    fams=[f for f in fam_n if keep(f)]
    ann=sum(1 for f in fams if any(m in annot for m in fam_rep.get(f,[])))
    n=len(fams)
    print(f"\n{title}: {n:,} families")
    print(f"   annotated:   {ann:>10,}  {100.0*ann/n:6.2f}%")
    print(f"   unannotated: {n-ann:>10,}  {100.0*(n-ann)/n:6.2f}%")
    return (title,n,ann)

print("\n=== 50% FAMILY TIER, dark fraction ===")
rows=[report("all families", lambda f: True),
      report(">=2 members", lambda f: fam_n[f]>=2),
      report(">=3 members and >=2 metagenomes", lambda f: fam_n[f]>=3 and len(fam_s[f])>=2),
      report(">=3 members and >=3 metagenomes", lambda f: fam_n[f]>=3 and len(fam_s[f])>=3)]

with open(W+"/results/family_dark_fraction.tsv","w") as fh:
    fh.write("filter\tfamilies\tannotated\tpct_annotated\tpct_unannotated\n")
    for t,n,a in rows:
        fh.write(f"{t}\t{n}\t{a}\t{100.0*a/n:.2f}\t{100.0*(n-a)/n:.2f}\n")
print("\nwrote", W+"/results/family_dark_fraction.tsv")
PYEOF
