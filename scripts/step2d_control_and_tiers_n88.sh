#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 8G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step2d_n88.log
# Step 2d n88: no new search. The new-cohort query set is catalog/subsample_200k.faa
# with UHM586.41010 queries removed; hits are per query against the fixed v1
# WF22 catalog, so the existing hit files restricted to the remaining queries
# equal a rerun. The WF22 positive control has no UHM586.41010 queries.
# Criteria as step2c/2d: --min-seq-id 0.95 (gene) or 0.5 (family), -c 0.8, --cov-mode 1.
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 << 'PYEOF'
import statistics
W="/bigdata/stajichlab/lshad003/wf_protein_space"
DROP="UHM586.41010__"
L={}; name=None
for line in open(W+"/catalog/subsample_200k.faa"):
    if line[0]==">": name=line[1:].strip(); L[name]=0
    else: L[name]+=len(line.strip())
q88=[q for q in L if not q.startswith(DROP)]
print("new-cohort queries: 89 %d, n88 %d (UHM586.41010 removed %d)"%(len(L),len(q88),len(L)-len(q88)))
def hits(p): return set(l.split("\t",1)[0] for l in open(p))
ctrlq=sum(1 for l in open(W+"/catalog/wf22_control_200k.faa") if l[0]==">")
ch=hits(W+"/results/step2d_control_hits.m8")
print("CONTROL (WF22 vs 95%% catalog, cov-mode 1): %d/%d = %.2f%%  [unchanged by n88]"%(len(ch),ctrlq,100.0*len(ch)/ctrlq))
s88=set(q88)
h50=hits(W+"/results/step2d_tier50_hits.m8")&s88
h95=hits(W+"/results/step2c_hits.m8")&s88
print("NEW vs 50%% tier (family level, cov-mode 1): %d/%d = %.2f%%"%(len(h50),len(q88),100.0*len(h50)/len(q88)))
print("NEW vs 95%% tier (gene level, cov-mode 1): %d/%d = %.2f%%"%(len(h95),len(q88),100.0*len(h95)/len(q88)))
m=[L[q] for q in q88 if q in h95]; u=[L[q] for q in q88 if q not in h95]
print("mapped:   n=%d median_aa=%d mean_aa=%d"%(len(m),statistics.median(m),sum(m)/len(m)))
print("unmapped: n=%d median_aa=%d mean_aa=%d"%(len(u),statistics.median(u),sum(u)/len(u)))
PYEOF
echo "DONE step2d_n88"
