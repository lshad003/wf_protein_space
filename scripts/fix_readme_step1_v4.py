import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 1\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 2\b.*$", t, re.M))
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
shutil.copy(f, f.replace("README.md", "README.pre_step1_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step1_v4.md").read() + t[b[0].start():]); print("written")
