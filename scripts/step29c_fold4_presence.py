import gzip, csv
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
NAME = {"UHM125.35769__k141_151130_2": "DivIVA-like", "UHM636.41030__k141_332084_5": "Cse2/CasB-like",
        "UHM145.35772__k141_106404_54": "HTH-like (weak)", "UHM145.35771__k141_106979_4": "helical repeat-like (weak)"}
md = {r["stem"]: r for r in csv.DictReader(open(W + "/metadata/samples_89_treatment.tsv"), delimiter="\t")}
mem = {q: [] for q in NAME}
with open(W + "/results/unk_direct30_clusters.tsv") as fh:
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if len(f) != 2: raise SystemExit("unexpected cluster line: %s" % f)
        if f[0] in mem: mem[f[0]].append(f[1])
sup = set(open(W + "/results/supported_reps_95.txt").read().split())
want = {m: q for q in NAME for m in mem[q] if m in sup}
cnt, S = {}, []
with gzip.open(W + "/results/count_matrix_primary.tsv.gz", "rt") as fh:
    S = fh.readline().rstrip("\n").split("\t")[1:]
    for l in fh:
        k = l[:l.find("\t")]
        if k in want:
            cnt[k] = [int(float(x)) for x in l.rstrip("\n").split("\t")[1:]]
            if len(cnt) == len(want): break
out = open(W + "/results/step29_fold4_presence.tsv", "w")
out.write("family_rep\tlabel\tsample\tanimal\tcohort\tegg_mass_cohort\ttreatment\tassembled\treads_rep95\n")
for q in NAME:
    asm = {m.split("__")[0] for m in mem[q]}
    r95 = [m for m in mem[q] if m in want]
    rd = [sum(cnt[m][i] for m in r95) for i in range(len(S))] if r95 else [0] * len(S)
    rdpos = {S[i] for i in range(len(S)) if rd[i] >= 5}
    print("\n== %s  %s | %d member proteins | assembled in %d samples, %d animals | supported rep95: %s; reads>=5 in %d samples, %d animals" % (
        q, NAME[q], len(mem[q]), len(asm), len({s.split(".")[0] for s in asm}), r95 or "none", len(rdpos), len({s.split(".")[0] for s in rdpos})))
    for key, lab in (("cohort", "cohort"), ("egg_mass_cohort", "egg mass"), ("treatment", "treatment")):
        t = {}
        for s in S:
            a = t.setdefault(md.get(s, {}).get(key, "NA"), [0, 0, 0]); a[2] += 1; a[0] += s in asm; a[1] += s in rdpos
        print("  %-9s %s" % (lab, " | ".join("%s asm %d, reads %d /%d" % (k, v[0], v[1], v[2]) for k, v in sorted(t.items()))))
    for i, s in enumerate(S):
        m = md.get(s, {})
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n" % (q, NAME[q], s, s.split(".")[0], m.get("cohort", "NA"), m.get("egg_mass_cohort", "NA"), m.get("treatment", "NA"), s in asm, rd[i]))
out.close()
