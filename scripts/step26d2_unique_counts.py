# Step 26d2: unique-sequence sizes of the 30% unknown families (step26a).
# For each 30% family in results/unk_clusters_30.tsv, take the members of its 95%
# clusters (results/clusters_95_n88.tsv) that are complete ORFs (partial=00) and
# >= 35 aa, and count unique sequences (exact duplicates removed, terminal * stripped).
# Sequences and flags are read from the per-sample Prodigal aa.fa.gz, the same source
# as Step 9 and step26a2. catalog/input_pep holds only the 45 WF23/WF24 stems, not
# WF22; for those 45 stems every sequence used here is checked identical to input_pep.
# Represented in NMPFams2 (metag) = any 95% rep of the family has a hit in
# results/nmpfams2_metag_rep_hits.tsv (E <= 1e-3, >=30% id, -c 0.7 --cov-mode 0).
# Read-only on every _n88 file.
import gzip, os, hashlib, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OLD="/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"
NEW=W+"/catalog/by_assembly_new"
PEP=W+"/catalog/input_pep"
OUT=W+"/results/unk_family_sizes_30_unique.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")

fam={}                                  # 95% rep -> 30% family rep
for line in open(W+"/results/unk_clusters_30.tsv"):
    f,r=line.rstrip("\n").split("\t"); fam[r]=f
mem={}                                  # member -> 95% rep
for line in open(W+"/results/clusters_95_n88.tsv"):
    r,m=line.rstrip("\n").split("\t")
    if r in fam: mem[m]=r
print("families:", len(set(fam.values())), "| member proteins:", len(mem), flush=True)

def fasta(fh):
    h=None; s=[]
    for line in fh:
        if line[0]==">":
            if h is not None: yield h,"".join(s)
            h=line[1:].rstrip("\n"); s=[]
        else: s.append(line.strip())
    if h is not None: yield h,"".join(s)

stems=sorted({m.split("__",1)[0] for m in mem})
good={}                                 # member -> md5 of sequence, complete >=35 aa only
seen=set(); checked=0; mismatch=0
for i,s in enumerate(stems,1):
    p = f"{NEW}/{s}.aa.fa.gz" if os.path.exists(f"{NEW}/{s}.aa.fa.gz") else f"{OLD}/{s}.aa.fa.gz"
    loc={}
    with gzip.open(p,"rt") as fh:
        for h,seq in fasta(fh):
            f=h.split(" # ")
            pid=s+"__"+f[0].split()[0]
            if pid not in mem: continue
            seen.add(pid)
            j=h.find("partial=")
            code=h[j+8:j+10] if j>=0 else "NA"
            seq=seq.rstrip("*")
            if code=="00" and len(seq)>=35:
                good[pid]=hashlib.md5(seq.encode()).digest(); loc[pid]=seq
    if os.path.exists(f"{PEP}/{s}.fasta"):
        with open(f"{PEP}/{s}.fasta") as fh:
            for h,seq in fasta(fh):
                k=h.split()[0]
                if k in loc:
                    checked+=1
                    if seq.rstrip("*")!=loc[k]: mismatch+=1
    if i%10==0: print(" ",i,"stems done", flush=True)
miss=[m for m in mem if m not in seen]
print("members without a Prodigal record:", len(miss), flush=True)
if miss: raise SystemExit("ABORT: missing members")
print(f"input_pep cross-check: {checked:,} sequences compared, {mismatch} differ", flush=True)
if mismatch: raise SystemExit("ABORT: sequence mismatch against input_pep")

ncomp=collections.Counter(); uniq=collections.defaultdict(set)
for m,d in good.items():
    f=fam[mem[m]]; ncomp[f]+=1; uniq[f].add(d)
nuniq={f:len(v) for f,v in uniq.items()}
nrep=collections.Counter(fam.values())
hit={l.split("\t",1)[0] for l in open(W+"/results/nmpfams2_metag_rep_hits.tsv")}
famhit={f for r,f in fam.items() if r in hit}
tc=sum(ncomp.values()); tu=sum(nuniq.values())
print(f"complete >=35 aa proteins {tc:,}; unique within families {tu:,} ({100.0*tu/tc:.2f}%)")

with open(OUT,"w") as fh:
    fh.write("family_rep\tn_reps_95\tn_proteins_complete35\tn_unique_complete35\tnmpfams2_metag\n")
    for f in sorted(nrep, key=lambda x:(-nuniq.get(x,0),x)):
        fh.write(f"{f}\t{nrep[f]}\t{ncomp[f]}\t{nuniq.get(f,0)}\t{int(f in famhit)}\n")
print("wrote", OUT)

print(f"\n{'count (complete >=35 aa)':28} {'>=25':>9} {'>=100':>9}")
for lab,c in (("proteins (step26a2)",ncomp),("unique sequences",nuniq)):
    print(f"{lab:28} " + " ".join(f"{sum(1 for f in nrep if c.get(f,0)>=t):>9,}" for t in (25,100)))
for t in (25,100):
    fs=[f for f in nrep if nuniq.get(f,0)>=t]; r=sum(1 for f in fs if f in famhit)
    print(f">= {t} unique: {len(fs):,} families, represented in NMPFams2 metag {r:,} ({100.0*r/len(fs):.2f}%)" if fs else f">= {t} unique: 0")
print("DONE step26d2")
