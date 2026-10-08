# Step 16u n88: both count matrices at 88 metagenomes, without remapping.
# Drops the UHM586.41010 column and keeps rows in results/supported_reps_95_n88.txt.
# The primary matrix already exists (step25a); it is re-derived here and must match.
# Column sums are checked against results/count_matrix_{primary,mapq10}.colsums.tsv.
import gzip, sys
W="/bigdata/stajichlab/lshad003/wf_protein_space/results/"
DROP="UHM586.41010"
sup=set(open(W+"supported_reps_95_n88.txt").read().split())
for tag in ("primary","mapq10"):
    ref={l.split("\t")[0]:int(l.split("\t")[1]) for l in open(W+"count_matrix_%s.colsums.tsv"%tag) if not l.startswith("stem")}
    fi=gzip.open(W+"count_matrix_%s.tsv.gz"%tag,"rt")
    h=fi.readline().rstrip("\n").split("\t")
    k=[i for i,x in enumerate(h) if x!=DROP]
    out=gzip.open(W+"count_matrix_%s_n88.check.tsv.gz"%tag if tag=="primary" else W+"count_matrix_mapq10_n88.tsv.gz","wt")
    out.write("\t".join(h[i] for i in k)+"\n")
    tot=[0]*len(h); nz89=0; nz88=0; rows=0
    for ln in fi:
        f=ln.rstrip("\n").split("\t"); v=[int(x) for x in f[1:]]
        for j,x in enumerate(v): tot[j+1]+=x
        if sum(v)>0: nz89+=1
        if f[0] in sup:
            rows+=1
            vv=[f[i] for i in k[1:]]
            if any(x!="0" for x in vv): nz88+=1
            out.write(f[0]+"\t"+"\t".join(vv)+"\n")
    out.close()
    bad=sum(1 for j,s in enumerate(h) if j>0 and tot[j]!=ref[s])
    print("%s: 89 genes with any count %d | n88 rows %d, genes with any count %d | column sums vs colsums.tsv mismatches %d"
          %(tag,nz89,rows,nz88,bad))
    if bad or rows!=2067011: sys.exit("ABORT %s"%tag)
print("DONE step16u_n88")
