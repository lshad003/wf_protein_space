import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 10\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 11\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["6,171,602", "48.29%", "36.83%", "90.53%", "65.1%", "40.4%", "76.03%", "55.40%", "1,020,958", "565,611"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 10")
shutil.copy(f, f.replace("README.md", "README.pre_step10_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step10_v4.md").read() + t[b[0].start():]); print("written")
