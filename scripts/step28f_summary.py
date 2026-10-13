# Step 28f summary: % of queries with any Foldseek hit at E <= 1e-1, 1e-2, 1e-3, 1e-5, 1e-10,
# for EU740 (results/prostt5_pilot, step28e/28e2), length-matched K740 and reversed EU740 decoys
# (results/prostt5_calib, step28f), separately vs PDB and CATH50. Denominator: sequences in each query FASTA.
# Output: results/prostt5_calib_summary.tsv; refuses to overwrite.
import os
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/prostt5_calib_summary.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")
SETS=(("EU740",W+"/results/prostt5_pilot_eu740.faa",W+"/results/prostt5_pilot/eu740_vs_%s.tsv"),
      ("K740_lenmatched",W+"/results/prostt5_calib_known_lenmatched.faa",W+"/results/prostt5_calib/known740_vs_%s.tsv"),
      ("EU740_reversed",W+"/results/prostt5_decoy_eu740rev.faa",W+"/results/prostt5_calib/decoy740_vs_%s.tsv"))
TH=(1e-1,1e-2,1e-3,1e-5,1e-10)
rows=[("set","db","n_queries")+tuple("pct_E<=%g"%t for t in TH)+tuple("n_E<=%g"%t for t in TH)]
for lab,faa,tpl in SETS:
    q={l[1:].split()[0] for l in open(faa) if l[0]==">"}
    for db in ("PDB","CATH50"):
        best={}
        for l in open(tpl%db):
            f=l.split("\t")
            if f[0] not in q: raise SystemExit("ABORT: hit query not in "+faa)
            e=float(f[2]); best[f[0]]=min(e,best.get(f[0],e))
        n=[sum(1 for e in best.values() if e<=t) for t in TH]
        rows.append((lab,db,len(q))+tuple("%.2f"%(100.0*k/len(q)) for k in n)+tuple(n))
with open(OUT,"w") as fh:
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
for r in rows: print("\t".join(map(str,r)))
