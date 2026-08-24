# Step 16h: extract CDS of the 2,069,453 supported 95% representatives
# from the 89 renamed per-stem CDS files into one mapping reference.
import glob
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
ids = set(open(W + "/results/supported_reps_95.txt").read().split())
assert len(ids) == 2069453, len(ids)
files = sorted(glob.glob(W + "/catalog/input_cds/*.fasta"))
print("input files:", len(files), flush=True)
out = open(W + "/catalog/db/LsPS_CDS_95_supported.fasta", "w")
seen = set()
keep = False
for i, f in enumerate(files, 1):
    for line in open(f):
        if line[0] == ">":
            name = line[1:].split()[0]
            keep = name in ids
            if keep:
                seen.add(name)
        if keep:
            out.write(line)
    if i % 10 == 0:
        print("files:", i, "found:", len(seen), flush=True)
out.close()
miss = ids - seen
lines = ["reference: supported reps only, >=3 members >=2 metagenomes, 95%% tier",
         "expected\t2069453", "found\t%d" % len(seen), "missing\t%d" % len(miss)]
with open(W + "/results/supported_cds_extract.summary.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
if miss:
    print("first missing:", sorted(miss)[:5])
print("EXTRACT:", "PASS" if not miss else "FAIL, do not index")
