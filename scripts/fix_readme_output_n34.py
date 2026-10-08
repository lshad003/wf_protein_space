f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
E = [
("results/limma_WF23_*.csv, results/limma_WF24_*.csv,",
 "results/limma_WF23_*.csv, results/limma_WF24_n34_*.csv (results/limma_WF24_treatment*.csv\n"
 "are the superseded 36-animal run),"),
("metadata/wf22_design.tsv, metadata/wf24_treatment_key.tsv",
 "metadata/wf22_design.tsv, metadata/wf24_treatment_key.tsv, metadata/wf24_excluded.tsv"),
]
for o, n in E: print(t.count(o), "|", o[:60])
if any(t.count(o) != 1 for o, n in E): raise SystemExit("ABORT: nothing written")
for o, n in E: t = t.replace(o, n)
open(f, "w").write(t); print("written")
