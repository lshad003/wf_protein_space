import shutil
W = "/bigdata/stajichlab/lshad003/wf_protein_space/"
t = open(W + "README.md").read()
s, e = "**WF22.** Three egg masses", "| File | Purpose |"
if t.count(s) != 1: raise SystemExit("ABORT: WF22 paragraph not found once")
i = t.index(s); j = t.index(e, i)
old = t[i:j]
if "**Exclusions.**" not in old or "**WF24.**" not in old: raise SystemExit("ABORT: block not as expected")
shutil.copy(W + "README.md", W + "README.pre_cohorts_v3.md")
open(W + "README.md", "w").write(t[:i] + open(W + "scripts/readme_cohorts_v3.md").read() + t[j:])
print("written")
