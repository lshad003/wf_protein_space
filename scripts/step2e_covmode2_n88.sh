#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 8G --time 1:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step2e_n88.log
# Step 2e n88: no new search; as step2d_n88, the existing step2e hit files are
# restricted to new-cohort queries outside UHM586.41010. Controls are WF22 only.
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 << 'PYEOF'
W="/bigdata/stajichlab/lshad003/wf_protein_space"
DROP="UHM586.41010__"
newq=[l[1:].strip() for l in open(W+"/catalog/subsample_200k.faa") if l[0]==">"]
s88=set(q for q in newq if not q.startswith(DROP))
ctrlq=sum(1 for l in open(W+"/catalog/wf22_control_200k.faa") if l[0]==">")
def hits(n): return set(l.split("\t",1)[0] for l in open(W+"/results/step2e_%s.m8"%n))
for mode,runs in (("cov-mode 2 (coverage of QUERY)",("control_cov2","new95_cov2","new50_cov2")),
                  ("cov-mode 0 (both sequences, strictest)",("control_cov0","new95_cov0"))):
    print("==== %s ===="%mode)
    for r in runs:
        h=hits(r)
        if r.startswith("control"): Q=ctrlq
        else: h&=s88; Q=len(s88)
        print("%-28s %7d/%7d = %6.2f%%"%(r,len(h),Q,100.0*len(h)/Q))
    print()
PYEOF
echo "DONE step2e_n88"
