import re, csv, collections, shutil
W = "/bigdata/stajichlab/lshad003/wf_protein_space/"
em = collections.Counter(r["egg_mass"] for r in csv.DictReader(open(W + "metadata/samples_88_treatment.tsv"), delimiter="\t") if r["cohort"] == "WF24")
print("WF24 egg masses:", dict(em))
if sorted(em.values()) != [3, 32]: raise SystemExit("ABORT: egg mass split is not 32/3, nothing written")
t = open(W + "README.md").read()
if t.count("## Sample set") != 1: raise SystemExit("ABORT: Sample set heading not found once")
m = list(re.finditer(r"^#+ Step 1\b.*$", t, re.M))
print("Step 1 headings:", [x.group(0) for x in m])
if len(m) != 1: raise SystemExit("ABORT: Step 1 heading not found exactly once")
new = open(W + "scripts/readme_header_v2.md").read()
shutil.copy(W + "README.md", W + "README.pre_header_v2.md")
open(W + "README.md", "w").write(new + t[m[0].start():])
print("written; old copy in README.pre_header_v2.md")
