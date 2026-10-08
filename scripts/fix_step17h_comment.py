f = "/bigdata/stajichlab/lshad003/wf_protein_space/scripts/step17h_limma_wf2324.R"
t = open(f).read()
old, new = "1 UHM520.7724,", "1 UHM520.7734,"
c = t.count(old)
print("found", c, "expected 1")
if c != 1: raise SystemExit("ABORT: nothing written")
open(f, "w").write(t.replace(old, new)); print("written")
