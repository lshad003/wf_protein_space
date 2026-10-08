f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
old, new = "| WF22 | 44 | 19 |", "| WF22 | 44 | 15 |"
c = t.count(old)
print("found", c, "expected 1")
if c != 1: raise SystemExit("ABORT: nothing written")
open(f, "w").write(t.replace(old, new)); print("written")
