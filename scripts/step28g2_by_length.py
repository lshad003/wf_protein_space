# Step 28g2: Foldseek hit rate by query length bin. Hit = best E over PDB and CATH50 <= 1e-3
# (main threshold; reversed decoys 0% in step28f). Length = residues excluding '*' stop characters.
# Bins: <100, 100-199, 200-399, >=400 aa, plus all. Sets: EU15580 (step28g), EU740 pilot (step28e/28e2),
# K740_lenmatched and EU740_reversed (step28f). Extra rows: EU740 pilot families found in EU15580
# (by representative ID), with their hit rate from the EU15580 search. Denominator: sequences in each FASTA.
# Output: results/prostt5_by_length.tsv; refuses to overwrite.
import os
W="/bigdata/stajichlab/lshad003/wf_protein_space"
OUT=W+"/results/prostt5_by_length.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: "+OUT+" exists")
EMAX=1e-3
SETS=(("EU15580",W+"/results/prostt5_main_eu15580.faa",W+"/results/prostt5_main/eu15580_vs_%s.tsv",15580),
      ("EU740",W+"/results/prostt5_pilot_eu740.faa",W+"/results/prostt5_pilot/eu740_vs_%s.tsv",740),
      ("K740_lenmatched",W+"/results/prostt5_calib_known_lenmatched.faa",W+"/results/prostt5_calib/known740_vs_%s.tsv",740),
      ("EU740_reversed",W+"/results/prostt5_decoy_eu740rev.faa",W+"/results/prostt5_calib/decoy740_vs_%s.tsv",740))
BINS=(("<100",0,100),("100-199",100,200),("200-399",200,400),(">=400",400,10**9))
def lengths(faa):
    L={}; h=None
    for l in open(faa):
        if l[0]==">": h=l[1:].split()[0]; L[h]=0
        else: L[h]+=len(l.strip().replace("*",""))
    return L
def hits(tpl,q):
    best={}
    for db in ("PDB","CATH50"):
        for l in open(tpl%db):
            f=l.split("\t")
            if f[0] not in q: raise SystemExit("ABORT: hit query not in query set: "+tpl%db)
            e=float(f[2]); best[f[0]]=min(e,best.get(f[0],e))
    return {x for x,e in best.items() if e<=EMAX}
def binrows(lab,L,H):
    out=[]
    for b,lo,hi in BINS+(("all",0,10**9),):
        ids=[x for x in L if lo<=L[x]<hi]; k=sum(1 for x in ids if x in H)
        out.append((lab,b,len(ids),k,"%.2f"%(100.0*k/len(ids)) if ids else "NA"))
    return out
rows=[("set","len_bin_aa","n_queries","n_hit_E<=1e-3_PDB_or_CATH50","pct_hit")]
D={}
for lab,faa,tpl,n in SETS:
    L=lengths(faa)
    if len(L)!=n: raise SystemExit("ABORT: %s has %d sequences, expected %d"%(faa,len(L),n))
    H=hits(tpl,L); D[lab]=(L,H)
    rows+=binrows(lab,L,H)
L15,H15=D["EU15580"]; L740,H740=D["EU740"]
ov=set(L740)&set(L15)
if any(L740[x]!=L15[x] for x in ov): raise SystemExit("ABORT: overlap IDs with different lengths")
rows+=binrows("EU740_in_EU15580_mainrun",{x:L15[x] for x in ov},H15)
rows+=binrows("EU740_in_EU15580_pilotrun",{x:L740[x] for x in ov},H740)
with open(OUT,"w") as fh:
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
for r in rows: print("\t".join(map(str,r)))
print("EU740 pilot families in EU15580: %d of %d"%(len(ov),len(L740)))
print("of these, hit in pilot run: %d; hit in main run: %d; hit in both: %d"%(len(ov&H740),len(ov&H15),len(ov&H740&H15)))
