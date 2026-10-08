import re, shutil, sys
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
start, end, md = sys.argv[1], sys.argv[2], sys.argv[3]
t = open(f).read()
a = list(re.finditer(r"^#+ Step " + start + r"\b.*$", t, re.M))
b = list(re.finditer(r"^#+ " + end + r".*$", t, re.M))
print("start:", [x.group(0) for x in a], "| end:", [x.group(0) for x in b])
if len(a) != 1 or len(b) != 1 or b[0].start() < a[0].start(): raise SystemExit("ABORT: headings, nothing written")
old = t[a[0].start():b[0].start()]
for n in sys.argv[4:]:
    if n not in old: raise SystemExit("ABORT: " + n + " not in current section")
shutil.copy(f, f.replace("README.md", "README.pre_step" + start + "_v4.md"))
open(f, "w").write(t[:a[0].start()] + open(md).read() + t[b[0].start():]); print("written")
