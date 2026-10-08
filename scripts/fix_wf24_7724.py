import shutil
M = "/bigdata/stajichlab/lshad003/wf_protein_space/metadata/"
edits = [
 ("samples_89_for_jason.tsv", "UHM520.7724", "UHM520.7734", 5),
 ("wf24_treatment_key.tsv", "UHM520.7724", "UHM520.7734", 1),
 ("wf24_treatment_key.tsv", "check vs WF23 UHM520.7734 possible typo",
  "typo corrected, confirmed by J. Dallas 2026-10-01", 1)]
text = {f: open(M + f).read() for f in {e[0] for e in edits}}
ok = True
for f, old, new, n in edits:
    c = text[f].count(old)
    print(f, repr(old), "found", c, "expected", n)
    if c != n: ok = False
    else: text[f] = text[f].replace(old, new)
if not ok:
    print("ABORT: counts differ, nothing written"); raise SystemExit(1)
for f, t in text.items():
    shutil.copy(M + f, M + f + ".bak_20261006")
    open(M + f, "w").write(t)
    print("written", f)
