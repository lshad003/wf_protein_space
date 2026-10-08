f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
o = "The catalog includes only metagenomes in the final experimental design.\n"
print("found", t.count(o), "expected 1")
if t.count(o) != 1: raise SystemExit("ABORT: nothing written")
open(f, "w").write(t.replace(o, "")); print("written")
