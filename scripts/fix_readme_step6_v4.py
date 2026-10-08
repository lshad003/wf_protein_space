import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 6\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 7\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["422,857", "431,504", "756,692", "1.75-fold", "1.24-fold"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 6")
shutil.copy(f, f.replace("README.md", "README.pre_step6_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step6_v4.md").read() + t[b[0].start():]); print("written")
