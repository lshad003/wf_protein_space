# Step 26a2: protein-level sizes of the 30% unknown families (step26a).
# For each 30% family in results/unk_clusters_30.tsv, sum the members of its
# constituent 95% clusters from results/clusters_95_n88.tsv (n88: no UHM586.41010
# members). Two counts per family:
#   all       every member protein
#   complete  members with Prodigal partial=00 and >= 35 aa
# Flags are read from the per-sample aa.fa.gz as in Step 9. Length of a complete
# ORF is (end - start + 1)/3 - 1 (stop codon not counted), same as len(seq.rstrip("*")).
# Writes results/unk_family_sizes_30_raw.tsv (per family) and prints thresholds.
# Read-only on every _n88 file.
import gzip, os, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
OUT=W+"/results/unk_family_sizes_30_raw.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")

fam={}                                  # 95% rep -> 30% family rep
for line in open(W+"/results/unk_clusters_30.tsv"):
    f,r=line.rstrip("\n").split("\t"); fam[r]=f
print("95% reps in 30% families:", len(fam), "| families:", len(set(fam.values())), flush=True)

mem={}                                  # member -> 95% rep
for line in open(W+"/results/clusters_95_n88.tsv"):
    r,m=line.rstrip("\n").split("\t")
    if r in fam: mem[m]=r
reps_seen=len(set(mem.values()))
print("member proteins:", len(mem), "| 95% reps with >=1 member in n88:", reps_seen, flush=True)

stems=sorted({m.split("__",1)[0] for m in mem})
print("stems with members:", len(stems), flush=True)
good={}                                 # member -> True if partial=00 and >=35 aa
for i,s in enumerate(stems,1):
    p = f"{NEW}/{s}.aa.fa.gz" if os.path.exists(f"{NEW}/{s}.aa.fa.gz") else f"{OLD}/{s}.aa.fa.gz"
    with gzip.open(p,"rt") as fh:
        for line in fh:
            if line[0]!=">": continue
            f=line[1:].split(" # ")
            pid=s+"__"+f[0].split()[0]
            if pid not in mem: continue
            j=line.find("partial=")
            code=line[j+8:j+10] if j>=0 else "NA"
            aa=(int(f[2])-int(f[1])+1)//3-1
            good[pid]=(code=="00" and aa>=35)
    if i%10==0: print(" ",i,"stems done", flush=True)
miss=[m for m in mem if m not in good]
print("members without a Prodigal header:", len(miss), flush=True)
if miss: raise SystemExit("ABORT: missing flags")

nrep=collections.Counter(fam.values())
nall=collections.Counter(); ncomp=collections.Counter()
for m,r in mem.items():
    nall[fam[r]]+=1
    if good[m]: ncomp[fam[r]]+=1
tot_all=sum(nall.values()); tot_comp=sum(ncomp.values())
print(f"\nmember proteins: all {tot_all:,} | complete and >=35 aa {tot_comp:,} ({100.0*tot_comp/tot_all:.2f}%)")

with open(OUT,"w") as fh:
    fh.write("family_rep\tn_reps_95\tn_proteins_all\tn_proteins_complete35\n")
    for f in sorted(nrep, key=lambda x:(-nall[x],x)):
        fh.write(f"{f}\t{nrep[f]}\t{nall[f]}\t{ncomp[f]}\n")
print("wrote", OUT)

print(f"\nfamilies total: {len(nrep):,}")
print(f"{'count':28} {'>=3':>9} {'>=25':>9} {'>=100':>9}")
for lab,c in (("95% representatives",nrep),("proteins, all members",nall),("proteins, complete >=35 aa",ncomp)):
    print(f"{lab:28} " + " ".join(f"{sum(1 for f in nrep if c[f]>=t):>9,}" for t in (3,25,100)))
print("DONE step26a2")
