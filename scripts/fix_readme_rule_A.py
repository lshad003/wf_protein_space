f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
o = "Because the WF24 exclusion was made after clustering"
n = "The catalog includes only metagenomes in the final experimental design.\nBecause the WF24 exclusion was made after clustering"
print("found", t.count(o), "expected 1")
if t.count(o) != 1: raise SystemExit("ABORT: nothing written")
open(f, "w").write(t.replace(o, n)); print("written")
