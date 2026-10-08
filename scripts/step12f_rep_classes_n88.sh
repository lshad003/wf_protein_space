#!/bin/bash
#SBATCH -p stajichlab -c 2 --mem 20G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step12f_n88.log
# Step 12f n88: per-rep flags and K/KWP/U class over results/reps_95_n88.txt,
# supported flag from results/supported_reps_95_n88.txt. Class is a property of
# the representative, so every n88 class must equal its class in rep_classes.tsv.
# DUF criterion as step12f: Pfam family name starts DUF or UPF (results/pfam_nonduf_ids.txt).
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 << 'PYEOF'
import sys
from collections import Counter
W="/bigdata/stajichlab/lshad003/wf_protein_space"
def load(p): return set(open(p).read().split())
egg=load(W+"/results/eggnog_hit_ids_n88.txt")
pf=load(W+"/results/pfam_hit_ids_n88.txt")
nd=load(W+"/results/pfam_nonduf_ids.txt")
sup=load(W+"/results/supported_reps_95_n88.txt")
old={}; oldsup=Counter(); oldduf=0
with open(W+"/results/rep_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p=ln.rstrip("\n").split("\t"); old[p[0]]=p[5]
        if p[1]=="1":
            oldsup[p[5]]+=1
            if p[5]=="K" and p[4]=="0": oldduf+=1
print("89 supported K/KWP/U:",oldsup["K"],oldsup["KWP"],oldsup["U"]," DUF/UPF-only K:",oldduf," strict K:",oldsup["K"]-oldduf)
out=open(W+"/results/rep_classes_n88.tsv","w")
out.write("rep\tsupported\teggnog\tpfam\tpfam_nonduf\tclass\n")
tot=Counter(); stot=Counter(); n=0; mism=0; duf=0
for line in open(W+"/results/reps_95_n88.txt"):
    r=line.strip(); n+=1
    e=r in egg; p=r in pf
    c="K" if p else ("KWP" if e else "U")
    if old.get(r)!=c: mism+=1
    tot[c]+=1
    s=r in sup
    if s:
        stot[c]+=1
        if c=="K" and r not in nd: duf+=1
    out.write("%s\t%d\t%d\t%d\t%d\t%s\n"%(r,s,e,p,r in nd,c))
out.close()
ns=sum(stot.values())
print("reps n88:",n)
print("all K/KWP/U:",tot["K"],tot["KWP"],tot["U"])
print("sup K/KWP/U:",stot["K"],stot["KWP"],stot["U"]," supported total:",ns)
print("sup pct K/KWP/U: %.2f %.2f %.2f"%tuple(100.0*stot[k]/ns for k in ("K","KWP","U")))
print("sup DUF/UPF-only K:",duf," strict K:",stot["K"]-duf)
print("class differs from rep_classes.tsv:",mism)
ok = mism==0 and ns==len(sup)==2067011
print("CROSS-CHECK:","PASS" if ok else "MISMATCH, inspect before use")
if not ok: sys.exit(1)
PYEOF
[ $? = 0 ] && echo "DONE step12f_n88"
