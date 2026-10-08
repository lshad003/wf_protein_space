f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
o, n = "Abundance across all 89 metagenomes", "Abundance across all 88 metagenomes"
print("found", t.count(o), "expected 2")
if t.count(o) != 2: raise SystemExit("ABORT: nothing written")
open(f, "w").write(t.replace(o, n)); print("written")
