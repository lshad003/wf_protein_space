import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Conventions.*$", t, re.M))
if len(a) != 1: raise SystemExit("ABORT: Conventions heading not found once")
nxt = re.search(r"^#+ ", t[a[0].end():], re.M)
end = a[0].end() + nxt.start() if nxt else len(t)
print("removing", len(t[a[0].start():end].splitlines()), "lines")
shutil.copy(f, f.replace("README.md", "README.pre_conventions.md"))
open(f, "w").write(t[:a[0].start()].rstrip("\n") + "\n" + ("\n" + t[end:] if nxt else "")); print("written")
