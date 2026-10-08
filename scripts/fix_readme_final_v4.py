import re, shutil
f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = list(re.finditer(r"^#+ Planned steps.*$", t, re.M)); b = list(re.finditer(r"^#+ Software.*$", t, re.M))
if len(a) != 1 or len(b) != 1: raise SystemExit("ABORT: Planned steps or Software heading")
t = t[:a[0].start()] + t[b[0].start():]
for old, new in [(16, 15), (17, 16), (19, 17), (20, 18), (21, 19), (22, 20)]:
    t, n = re.subn(r"^(#+ Step )%d\b" % old, r"\g<1>%d" % new, t, flags=re.M)
    print("Step", old, "->", new, ":", n)
    if n != 1: raise SystemExit("ABORT: heading count, nothing written")
def flex(s): return r"\s+".join(re.escape(w) for w in s.split())
E = [("with Prodigal in metagenomic mode", "with Prodigal v2.6.3 in metagenomic mode"),
     ("MMseqs2 is pinned to 13-45111 and runs only on the epyc partition, since other nodes fail with an illegal instruction.", "")]
for o, n in E:
    t, c = re.subn(r"\s?" + flex(o) if n == "" else flex(o), n, t)
    print(c, "|", o[:50])
    if c != 1: raise SystemExit("ABORT: text not found once, nothing written")
shutil.copy(f, f.replace("README.md", "README.pre_final_v4.md"))
open(f, "w").write(t); print("written")
