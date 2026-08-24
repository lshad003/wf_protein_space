# Step 16j helper: id and length for every supported CDS, for normalization.
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
REF = W + "/catalog/db/LsPS_CDS_95_supported.fasta"
OUT = W + "/results/supported_cds_lengths.tsv"
name = None
L = 0
n = 0
out = open(OUT, "w")
for line in open(REF):
    if line[0] == ">":
        if name is not None:
            out.write("%s\t%d\n" % (name, L))
        name = line[1:].split()[0]
        L = 0
        n += 1
    else:
        L += len(line.strip())
if name is not None:
    out.write("%s\t%d\n" % (name, L))
out.close()
print("length table:", n, "sequences, expect 2069453")
assert n == 2069453, "sequence count mismatch, stop"
