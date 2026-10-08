import csv, collections
M = "/bigdata/stajichlab/lshad003/wf_protein_space/metadata/"
rows = [r for r in csv.DictReader(open(M + "samples_89_treatment.tsv"), delimiter="\t") if r["cohort"] == "WF22"]
print("WF22 rows:", len(rows))
by = collections.defaultdict(list)
for r in rows: by[r["stem"].split(".")[0]].append(r)
print("distinct stem prefixes:", len(by))
print("samples per prefix:", sorted(collections.Counter(len(v) for v in by.values()).items()))
for a, v in sorted(by.items()):
    print(a, len(v), "|", " ".join(x["stem"] for x in v), "| trt", sorted({x["treatment"] for x in v}),
          "| em", sorted({x["egg_mass_cohort"] for x in v}))
print("\n== wf22_design.tsv ==")
d = list(csv.DictReader(open(M + "wf22_design.tsv"), delimiter="\t"))
print("rows:", len(d), "columns:", list(d[0].keys()))
for col in d[0].keys():
    if any(k in col.lower() for k in ("animal", "frog", "id", "subject", "individual")):
        vals = collections.Counter(r[col] for r in d)
        print(col, "distinct:", len(vals), "| per value:", sorted(collections.Counter(vals.values()).items()))
