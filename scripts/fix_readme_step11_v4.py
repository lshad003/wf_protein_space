import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 11\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 12\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["2,980,052", "92.64%", "6.93%", "0.27%", "0.16%", "16.30%", "13.95%", "13.34%", "11.44%", "10.31%", "3.05%", "2.57%", "7.63%", "3.61%"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 11")
shutil.copy(f, f.replace("README.md", "README.pre_step11_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step11_v4.md").read() + t[b[0].start():]); print("written")
