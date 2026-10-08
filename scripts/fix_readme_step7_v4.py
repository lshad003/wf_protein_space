import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 7\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 8\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["1,169,677", "1,159,262", "1,295,364", "1,687,560", "30.3%", "11.7%", "6,171,602", "40,257,962", "436,425"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 7")
shutil.copy(f, f.replace("README.md", "README.pre_step7_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step7_v4.md").read() + t[b[0].start():]); print("written")
