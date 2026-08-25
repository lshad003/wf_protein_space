# Step 21b: gene neighbourhood of unknown families.
# For each unknown family with >=3 members, find the genes immediately
# before and after each member on its contig, and count which Pfam
# domains those neighbours carry. Conserved neighbours suggest function.
import collections, gzip
W="/bigdata/stajichlab/lshad003/wf_protein_space"

# family membership for unknown genes
fam={}
for ln in open(W+"/results/unk_clusters_50.tsv"):
    r,m=ln.rstrip("\n").split("\t")[:2]
    fam[m]=r
size=collections.Counter(fam.values())
keep={r for r,n in size.items() if n>=3}
print("unknown families >=3 members:",len(keep),flush=True)

# map every catalog member to its 95% representative
rep={}
for ln in open(W+"/results/clusters_95.tsv"):
    r,m=ln.rstrip("\n").split("\t")[:2]
    rep[m]=r
print("catalog members:",len(rep),flush=True)

# Pfam domain per representative, first domain per protein
dom={}
import glob
for f in glob.glob(W+"/results/pfam_hs/*.tblout"):
    for ln in open(f):
        if ln[0]=="#": continue
        p=ln.split()
        if p[0] not in dom: dom[p[0]]=p[2]
print("reps with a Pfam domain:",len(dom),flush=True)

def nb(g,off):
    stem,rest=g.split("__",1)
    i=rest.rfind("_")
    try: n=int(rest[i+1:])
    except ValueError: return None
    return "%s__%s_%d"%(stem,rest[:i],n+off)

counts=collections.defaultdict(collections.Counter)
ncontig=collections.Counter()
for m,f in fam.items():
    if f not in keep: continue
    ncontig[f]+=1
    seen=set()
    for off in (-1,1):
        x=nb(m,off)
        if x is None: continue
        r=rep.get(x)
        if r is None: continue
        d=dom.get(r)
        if d: seen.add(d)
    for d in seen: counts[f][d]+=1

out=open(W+"/results/unk_family_neighbours.tsv","w")
out.write("family\tn_members\ttop_domain\tn_with_domain\tpct_of_members\tsecond_domain\tn_second\n")
strong=0
rows=[]
for f in keep:
    c=counts[f]; n=ncontig[f]
    if not c:
        out.write("%s\t%d\tnone\t0\t0.00\t\t0\n"%(f,n)); continue
    top=c.most_common(2)
    pct=100.0*top[0][1]/n
    if pct>=50: strong+=1
    s=top[1] if len(top)>1 else ("",0)
    out.write("%s\t%d\t%s\t%d\t%.2f\t%s\t%d\n"%(f,n,top[0][0],top[0][1],pct,s[0],s[1]))
    rows.append((pct,top[0][1],f,n,top[0][0]))
out.close()

print("\nfamilies with a domain in >=50%% of members: %d of %d"%(strong,len(keep)))
print("\n=== top 20 unknown families by neighbour conservation (>=5 members) ===")
print("pct    hits  members  family                              top neighbour domain")
for pct,h,f,n,d in sorted([r for r in rows if r[3]>=5], reverse=True)[:20]:
    print("%-6.1f %-5d %-8d %-35s %s"%(pct,h,n,f,d))
print("\nwrote results/unk_family_neighbours.tsv")
