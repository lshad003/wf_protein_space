#!/bin/bash
#SBATCH -p stajichlab -c 2 --mem 20G --time 4:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step12f.log
# Step 12f: per-rep flags (supported, eggnog, pfam, pfam_nonduf) and provisional
# class K/KWP/U. DUF criterion: Pfam family name starts DUF or UPF (tblout col 3).
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
R=$WPS/results
export LC_ALL=C
cat $R/pfam_hs/*.tblout | grep -v '^#' \
| awk '$3 !~ /^DUF/ && $3 !~ /^UPF/ {print $1}' \
| sort -u -S 6G -T $WPS/tmp > $R/pfam_nonduf_ids.txt
echo "nonduf reps: $(wc -l < $R/pfam_nonduf_ids.txt) of pfam $(wc -l < $R/pfam_hit_ids.txt)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 << 'PYEOF'
W="/bigdata/stajichlab/lshad003/wf_protein_space"
def load(p): return set(open(p).read().split())
egg=load(W+"/results/eggnog_hit_ids.txt")
pf=load(W+"/results/pfam_hit_ids.txt")
nd=load(W+"/results/pfam_nonduf_ids.txt")
sup=load(W+"/results/supported_reps_95.txt")
out=open(W+"/results/rep_classes.tsv","w")
out.write("rep\tsupported\teggnog\tpfam\tpfam_nonduf\tclass\n")
from collections import Counter
tot=Counter(); stot=Counter(); n=0
for line in open(W+"/catalog/db/LsPS_AA_95_rep.fasta"):
    if line[0]!=">": continue
    r=line[1:].split()[0]; n+=1
    e=r in egg; p=r in pf
    c="K" if p else ("KWP" if e else "U")
    tot[c]+=1
    s=r in sup
    if s: stot[c]+=1
    out.write("%s\t%d\t%d\t%d\t%d\t%s\n"%(r,s,e,p,r in nd,c))
out.close()
print("reps:",n,"expect 6182117")
print("all K/KWP/U:",tot["K"],tot["KWP"],tot["U"],"expect 2713156 514128 2954833")
print("sup K/KWP/U:",stot["K"],stot["KWP"],stot["U"],"expect 1344403 138835 586215")
ok=(n==6182117 and tot["K"]==2713156 and tot["KWP"]==514128 and tot["U"]==2954833
    and stot["K"]==1344403 and stot["KWP"]==138835 and stot["U"]==586215)
print("CROSS-CHECK:","PASS" if ok else "MISMATCH, inspect before use")
PYEOF
echo "DONE step12f"
