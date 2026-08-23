#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 4 --mem 120gb --time=8:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step13.%j.log
# Taxonomic composition of the annotated representatives, from the eggNOG
# max_annot_lvl and eggNOG_OGs fields. Asks what kingdom-level groups the
# annotated fraction is drawn from, and whether eukaryotic assignment tracks
# cluster size, since host and dietary sequence would concentrate in poorly
# supported clusters.
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
$PY << 'PYEOF'
import glob, collections
W="/bigdata/stajichlab/lshad003/wf_protein_space"

files=sorted(glob.glob(W+"/results/eggnog2/*.emapper.annotations"))
hdr=None
for line in open(files[0]):
    if line.startswith("#query"):
        hdr=line.lstrip("#").rstrip("\n").split("\t"); break
print("columns:", hdr)
iOG = hdr.index("eggNOG_OGs")
iLV = hdr.index("max_annot_lvl")

def kingdom(ogs):
    # the root-level OG carries the broad clade after the pipe
    for og in ogs.split(","):
        if "|" in og:
            tax = og.split("|",1)[1]
            if tax.startswith("2|") or tax=="Bacteria": return "Bacteria"
            if tax.startswith("2759") or "Eukaryota" in tax: return "Eukaryota"
            if tax.startswith("2157") or "Archaea" in tax: return "Archaea"
            if "Viruses" in tax: return "Viruses"
    return "other"

tax=collections.Counter(); lvl=collections.Counter(); prot_tax={}
for f in files:
    for line in open(f):
        if line.startswith("#"): continue
        p=line.rstrip("\n").split("\t")
        k=kingdom(p[iOG]); tax[k]+=1; lvl[p[iLV]]+=1
        prot_tax[p[0]]=k
tot=sum(tax.values())
print(f"\n=== broad clade of annotated representatives (n={tot:,}) ===")
for k,v in tax.most_common():
    print(f"  {k:12} {v:>10,}  {100.0*v/tot:6.2f}%")

print("\n=== most frequent max_annot_lvl (top 15) ===")
for k,v in lvl.most_common(15):
    print(f"  {k:40} {v:>10,}  {100.0*v/tot:6.2f}%")

print("\nloading cluster sizes ...", flush=True)
size=collections.Counter()
for line in open(W+"/results/clusters_95.tsv"):
    size[line.split("\t",1)[0]]+=1

print("\n=== clade against cluster size (annotated only) ===")
bins=[(1,1,"singleton"),(2,2,"2"),(3,5,"3-5"),(6,20,"6-20"),(21,10**9,"21+")]
print(f"{'bin':12}" + "".join(f"{k:>12}" for k in ["Bacteria","Eukaryota","Archaea","other"]))
rows=[]
for lo,hi,lab in bins:
    c=collections.Counter()
    for p,k in prot_tax.items():
        s=size.get(p,0)
        if lo<=s<=hi: c[k]+=1
    n=sum(c.values())
    print(f"{lab:12}" + "".join(f"{100.0*c[k]/n if n else 0:11.2f}%" for k in ["Bacteria","Eukaryota","Archaea","other"]))
    rows.append((lab,n,c))

with open(W+"/results/taxonomy_summary.tsv","w") as fh:
    fh.write("category\tn\tBacteria\tEukaryota\tArchaea\tother\n")
    fh.write(f"all_annotated\t{tot}\t{tax['Bacteria']}\t{tax['Eukaryota']}\t{tax['Archaea']}\t{tax['other']}\n")
    for lab,n,c in rows:
        fh.write(f"size_{lab}\t{n}\t{c['Bacteria']}\t{c['Eukaryota']}\t{c['Archaea']}\t{c['other']}\n")
print("\nwrote", W+"/results/taxonomy_summary.tsv")
PYEOF
