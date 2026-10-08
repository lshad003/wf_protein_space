f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
a = "| scripts/step17h_limma_wf2324.R | WF23 and WF24 treatment contrasts |"
b = "| scripts/step22c_jaccard.R | Presence/absence composition and gene richness |"
c = "Output: logs step22b, step22c and step22d;"
E = [
(a, "| scripts/step17h_limma_wf2324.R | WF23 treatment contrasts; its WF24 run (36 animals) is superseded by step17i |\n"
    "| scripts/step17i_limma_wf24_n34.R | WF24 treatment contrasts, 34 animals, exclusions in metadata/wf24_excluded.tsv |"),
(b, b + "\n| scripts/step22e_permanova_wf24_n34.R | WF24 abundance and presence/absence composition, 34 animals |"),
(c, "Output: logs step22b, step22c, step22d and step22e;"),
]
for o, n in E: print(t.count(o), "|", o[:60])
if any(t.count(o) != 1 for o, n in E): raise SystemExit("ABORT: nothing written")
for o, n in E: t = t.replace(o, n)
open(f, "w").write(t); print("written")
