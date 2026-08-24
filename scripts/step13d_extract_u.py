# Step 13d helper: fasta of supported class-U reps (no eggNOG, no Pfam).
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
ids = set()
with open(W + "/results/rep_classes.tsv") as fh:
    next(fh)
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if p[1] == "1" and p[5] == "U":
            ids.add(p[0])
assert len(ids) == 586215, len(ids)
out = open(W + "/results/supported_U_reps.faa", "w")
keep = False
n = 0
for line in open(W + "/catalog/db/LsPS_AA_95_rep.fasta"):
    if line[0] == ">":
        keep = line[1:].split()[0] in ids
        if keep:
            n += 1
    if keep:
        out.write(line)
out.close()
print("wrote", n, "sequences, expect 586215")
assert n == 586215, n
