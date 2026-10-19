# Step 34j B3a2: write the Merizo domain PDBs that step34j_b3a_merizo.sh (job 30050461) failed to write.
# Its inline writer crashed (KeyError 'filename': Merizo 41d12fb names the column 'input' and repeats the header row
# before every model), leaving results/step34j/domains empty although Merizo completed for all 1,156 models >= 200 aa.
# Same rule as that writer, on the existing results/step34j/merizo/merizo_results.tsv (Merizo is not rerun):
# result column, domains separated by ",", discontinuous segments by "_"; results/step34j/domains/<model>__D<k>.pdb
# holds the ATOM lines of domain k from results/step34j/esmfold_pdb/<model>.pdb. Rows with ndom 0 (no result column)
# write nothing. Refuses to run if domains/ already holds PDBs.
import os
D = "/bigdata/stajichlab/lshad003/wf_protein_space/results/step34j/"
if any(f.endswith(".pdb") for f in os.listdir(D + "domains")): raise SystemExit("REFUSING: domains/ not empty")
h = None; rows = []
for l in open(D + "merizo/merizo_results.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "input": h = f; continue
    if l.strip(): rows.append(dict(zip(h, f)))
names = [os.path.basename(r["input"])[:-4] for r in rows]
if len(set(names)) != len(names): raise SystemExit("ABORT: duplicate models in merizo_results.tsv")
print("models in merizo_results.tsv: %d" % len(rows))
nd = nm = 0; exp = 0
for r, name in zip(rows, names):
    exp += int(r["ndom"])
    res = r.get("result", "").strip()
    if not res or res in ("NA", "-"): continue
    doms = res.split(",")
    if len(doms) != int(r["ndom"]): raise SystemExit("ABORT: ndom disagrees with result for " + name)
    lines = [l for l in open(D + "esmfold_pdb/%s.pdb" % name) if l.startswith("ATOM")]
    nm += 1
    for k, d in enumerate(doms, 1):
        keep = set()
        for seg in d.split("_"):
            a, b = seg.split("-"); keep |= set(range(int(a), int(b) + 1))
        out = [l for l in lines if int(l[22:26]) in keep]
        if not out: raise SystemExit("ABORT: empty domain %s D%d" % (name, k))
        with open(D + "domains/%s__D%d.pdb" % (name, k), "w") as o:
            o.writelines(out); o.write("END\n")
        nd += 1
print("models with domains %d of %d (ndom 0: %d); domain PDBs written %d (sum of ndom %d)" % (nm, len(rows), len(rows) - nm, nd, exp))
if nd != exp: raise SystemExit("ABORT: written count differs from sum of ndom")
