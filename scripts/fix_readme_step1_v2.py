import re, shutil
W = "/bigdata/stajichlab/lshad003/wf_protein_space/"
t = open(W + "README.md").read()
a = list(re.finditer(r"^#+ Step 1\b.*$", t, re.M))
b = list(re.finditer(r"^#+ Step 2\b.*$", t, re.M))
print("Step 1:", [x.group(0) for x in a], "| Step 2:", [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1 or b[0].start() < a[0].start():
    raise SystemExit("ABORT: headings not found exactly once, nothing written")
old = t[a[0].start():b[0].start()]
if "40,550,595" not in old or "UHM585" not in old: raise SystemExit("ABORT: Step 1 text not as expected")
new = open(W + "scripts/readme_step1_v2.md").read()
h = a[0].group(0).split("Step")[0]
new = new.replace("### Step 1", h + "Step 1", 1)
shutil.copy(W + "README.md", W + "README.pre_step1_v2.md")
open(W + "README.md", "w").write(t[:a[0].start()] + new + t[b[0].start():])
print("written")
