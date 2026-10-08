import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
o = ", since it changes recovery by more than thirty percentage points"
p = r",\s+".join(r"\s+".join(re.escape(w) for w in part.split()) for part in o.split(",")[1:])
t, c = re.subn(r"," + r"\s+" + p, "", t)
print("coverage clause found", c)
if c != 1: raise SystemExit("ABORT: coverage clause not found once")
a = list(re.finditer(r"^#+ Interpretation rules.*$", t, re.M))
if len(a) != 1: raise SystemExit("ABORT: Interpretation rules heading not found once")
nxt = re.search(r"^#+ ", t[a[0].end():], re.M)
end = a[0].end() + nxt.start() if nxt else len(t)
print("removing", len(t[a[0].start():end].splitlines()), "lines")
shutil.copy(f, f.replace("README.md", "README.pre_rules.md"))
open(f, "w").write(t[:a[0].start()].rstrip("\n") + "\n" + ("\n" + t[end:] if nxt else "")); print("written")
