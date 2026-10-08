import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 13\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 14\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["585,347", "470,748,714", "2.1.24", "430,299", "20.82%", "155,048", "7.50%", "63.92%", "76.63%"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 13")
shutil.copy(f, f.replace("README.md", "README.pre_step13_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step13_v4.md").read() + t[b[0].start():]); print("written")
