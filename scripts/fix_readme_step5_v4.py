import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 5\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 6\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["15.27%", "9.52%", "25.28%", "24.33%", "765,122", "1,166,662", "2,080,538"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 5")
shutil.copy(f, f.replace("README.md", "README.pre_step5_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step5_v4.md").read() + t[b[0].start():]); print("written")
