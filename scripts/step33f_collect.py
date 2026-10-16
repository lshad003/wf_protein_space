# Step 33f collect: every barrnap 18S_rRNA hit >= 500 bp from results/step33f_barrnap/<stem>.rrna.fa (88 samples),
# headers "<stem>|<barrnap header>" (barrnap 0.9: 18S_rRNA::<contig>:<start>-<end>(<strand>)).
# Output (refuses to overwrite): results/step33f_18S_ge500.fa
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = R + "step33f_18S_ge500.fa"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
stems = [l.split("\t")[0] for i, l in enumerate(open(W + "/metadata/samples_88_treatment.tsv")) if i]
n = collections.Counter(); L = []
with open(OUT, "w") as o:
    for s in stems:
        p = R + "step33f_barrnap/%s.rrna.fa" % s
        if not os.path.exists(p): raise SystemExit("ABORT: missing " + p)
        h, seq = None, []
        for l in list(open(p)) + [">"]:
            if l[0] == ">":
                if h is not None:
                    q = "".join(seq); kind = h.split("::")[0]; n[kind] += 1
                    if kind == "18S_rRNA":
                        n["18S_all"] += 1
                        if len(q) >= 500: o.write(">%s|%s\n%s\n" % (s, h, q)); L.append(len(q)); n["18S_ge500"] += 1
                h, seq = l[1:].strip(), []
            else: seq.append(l.strip())
print("barrnap hits by type: %s" % {k: v for k, v in n.items() if k.endswith("rRNA")})
print("18S hits: %d | >= 500 bp: %d" % (n["18S_all"], n["18S_ge500"]))
if L: L.sort(); print("18S >= 500 bp length: min %d, median %d, max %d" % (L[0], L[len(L) // 2], L[-1]))
print("wrote", OUT)
