import gzip, sys
W="/bigdata/stajichlab/lshad003/wf_protein_space"
COL=int(sys.argv[1]); TAG=sys.argv[2]
stems=[l.split("\t")[0] for l in open(W+"/results/cohort_map.tsv")]
stems=[s.strip() for s in stems]
assert len(stems)==89
reps=[l.strip() for l in open(W+"/results/supported_reps_95.txt")]
idx={r:i for i,r in enumerate(reps)}
n=len(reps); assert n==2069453
out=gzip.open(W+"/results/count_matrix_%s.tsv.gz"%TAG,"wt")
out.write("gene\t"+"\t".join(stems)+"\n")
cols=[]
for k,s in enumerate(stems,1):
    v=[0]*n
    for ln in gzip.open(W+"/results/counts/%s.counts.tsv.gz"%s,"rt"):
        p=ln.rstrip("\n").split("\t")
        i=idx.get(p[0])
        if i is not None: v[i]=int(p[COL])
    cols.append(v)
    print("read",k,s,flush=True)
tot=[0]*89; nz=0
for i in range(n):
    row=[cols[j][i] for j in range(89)]
    s=sum(row)
    if s>0: nz+=1
    for j in range(89): tot[j]+=row[j]
    out.write(reps[i]+"\t"+"\t".join(map(str,row))+"\n")
out.close()
print("genes:",n,"genes with any count:",nz)
print("column sums written to summary")
with open(W+"/results/count_matrix_%s.colsums.tsv"%TAG,"w") as fh:
    fh.write("stem\ttotal_counts\n")
    for j,s in enumerate(stems): fh.write("%s\t%d\n"%(s,tot[j]))
