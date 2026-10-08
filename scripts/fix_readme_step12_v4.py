import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Step 12\b.*$", t, re.M)); b = list(re.finditer(r"^#+ Step 13\b.*$", t, re.M))
print([x.group(0) for x in a], [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: nothing written")
old = t[a[0].start():b[0].start()]
for n in ["3,222,677", "52.22%", "2,067,011", "585,347", "28.32%", "51.47%", "55.40%", "1,343,020", "76,515", "138,644"]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current Step 12")
shutil.copy(f, f.replace("README.md", "README.pre_step12_v4.md"))
open(f, "w").write(t[:a[0].start()] + open("scripts/readme_step12_v4.md").read() + t[b[0].start():]); print("written")
