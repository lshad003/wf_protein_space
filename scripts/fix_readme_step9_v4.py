import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 9\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 10\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["46.89%", "40,257,962", "11.58%", "32.98%", "61.16%", "81.99%", "86.11%", "2,067,011"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 9")
shutil.copy(f, f.replace("README.md", "README.pre_step9_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step9_v4.md").read() + t[b[0].start():]); print("written")
