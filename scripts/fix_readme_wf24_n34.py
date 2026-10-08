f = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(f).read()
E = [
("3 egg mass 2, across 7 treatment codes with 5 animals each. Tank is confounded",
 "3 egg mass 2, across 7 treatment codes; codes 1 to 6 have 5 animals and code 7\n"
 "has 6. Two of these animals were dropped from the experiment by J. Dallas and\n"
 "are excluded from all treatment tests: UHM586, a code 4 control with elevated\n"
 "*Basidiobolus* reads from two treatments, and UHM590 (code 7, reason not\n"
 "recorded). WF24 treatment tests therefore use 34 animals, 4 controls and 5 per\n"
 "strain. In 2024 treated animals received at least two inoculations within the\n"
 "month. Tank is confounded"),
("where 36 animals are each sampled once, treatment explains 14.9% (p = 0.898).",
 "where 34 animals are each sampled once, treatment R2 is 0.173 (p = 0.678), at\n"
 "the chance level expected for 6 of 33 degrees of freedom."),
("and in WF24 14.8% (p = 0.951).",
 "and in WF24 R2 0.170 (p = 0.770, 34 animals, chance level)."),
("and one gene for the sixth.",
 "and three genes for the sixth (code 3, UHM516.7697; 34 animals)."),
]
bad = [(o, t.count(o)) for o, n in E if t.count(o) != 1]
for o, n in E: print(t.count(o), "|", o[:60])
if bad: raise SystemExit("ABORT: nothing written")
for o, n in E: t = t.replace(o, n)
open(f, "w").write(t); print("written")
