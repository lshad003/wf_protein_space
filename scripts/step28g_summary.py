# Step 28g summary: % and n of the 15,580 queries (results/prostt5_main_eu15580.faa) with any Foldseek hit
# at E <= 1e-1, 1e-2, 1e-3, 1e-5, 1e-10 vs PDB, CATH50 and either (best E over both databases).
# Main threshold E <= 1e-3 (reversed decoys 0% in step28f). Output: results/prostt5_main_summary.tsv.
import os
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/prostt5_main_summary.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")
q={l[1:].split()[0] for l in open(W+"/results/prostt5_main_eu15580.faa") if l[0]==">"}
if len(q)!=15580: raise SystemExit("ABORT: query set is not 15,580")
TH=(1e-1,1e-2,1e-3,1e-5,1e-10)
best={}
for db in ("PDB","CATH50"):
    b={}
    for l in open(W+"/results/prostt5_main/eu15580_vs_%s.tsv"%db):
        f=l.split("\t")
        if f[0] not in q: raise SystemExit("ABORT: hit query not in query set")
        e=float(f[2]); b[f[0]]=min(e,b.get(f[0],e))
    best[db]=b
best["either"]={x:min(best["PDB"].get(x,9e9),best["CATH50"].get(x,9e9)) for x in set(best["PDB"])|set(best["CATH50"])}
rows=[("db","n_queries")+tuple("pct_E<=%g"%t for t in TH)+tuple("n_E<=%g"%t for t in TH)]
for db in ("PDB","CATH50","either"):
    n=[sum(1 for e in best[db].values() if e<=t) for t in TH]
    rows.append((db,len(q))+tuple("%.2f"%(100.0*k/len(q)) for k in n)+tuple(n))
with open(OUT,"w") as fh:
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
for r in rows: print("\t".join(map(str,r)))
