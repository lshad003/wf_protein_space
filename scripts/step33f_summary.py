# Step 33f summary: composition of the 18S hits >= 500 bp on Tiara-eukaryotic contigs (results/step33f_18S_ge500.fa)
# from their top SILVA SSU 138.1 NR99 hit (vsearch --usearch_global, results/step33f_18S_silva.tsv; first of tied
# top hits). Groups from the SILVA path: non-Eukaryota top hit; Metazoa: Amphibia / Insecta / other; Fungi;
# Ciliophora; Opalinata; Stramenopiles (other); Viridiplantae (SILVA Chloroplastida); other protists (all other
# Eukaryota); no hit (< --id). Per cohort (metadata/samples_88_treatment.tsv): number of 18S hits and number of
# samples with >= 1 hit in the group. Descriptive only: cohort is confounded with sequencing run.
# Outputs (refuse to overwrite): results/step33f_18S_classified.tsv, results/step33f_composition.tsv
import os, re, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = [R + "step33f_18S_classified.tsv", R + "step33f_composition.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
coh = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: f = l.split("\t"); coh[f[0]] = f[1]
q = {}
for l in open(R + "step33f_18S_ge500.fa"):
    if l[0] == ">": h = l[1:].strip(); q[h] = 0
    else: q[h] += len(l.strip())
top, ties = {}, collections.Counter()
for l in open(R + "step33f_18S_silva.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] not in q: raise SystemExit("ABORT: unknown query " + f[0])
    ties[f[0]] += 1
    if f[0] not in top: top[f[0]] = f
if set(top) != set(q): raise SystemExit("ABORT: %d queries missing from vsearch output" % len(set(q) - set(top)))
GROUPS = ["Metazoa: Amphibia", "Metazoa: Insecta", "Metazoa: other", "Fungi", "Ciliophora", "Opalinata",
          "Stramenopiles (other)", "Viridiplantae", "other protists", "non-eukaryotic top hit", "no hit"]
def group(target):
    if target == "*": return "no hit", "-"
    tax = target.split(" ", 1)[1] if " " in target else ""; p = tax.split(";")
    if p[0] != "Eukaryota": return "non-eukaryotic top hit", tax
    if "Metazoa" in p: return ("Metazoa: Amphibia" if "Amphibia" in p else "Metazoa: Insecta" if "Insecta" in p else "Metazoa: other"), tax
    if "Fungi" in p: return "Fungi", tax
    if "Ciliophora" in p: return "Ciliophora", tax
    if "Opalinata" in p: return "Opalinata", tax
    if "Stramenopiles" in p: return "Stramenopiles (other)", tax
    if "Chloroplastida" in p: return "Viridiplantae", tax
    return "other protists", tax
hits = collections.defaultdict(collections.Counter); samp = collections.defaultdict(lambda: collections.defaultdict(set))
ids = collections.defaultdict(list); lowest = collections.defaultdict(collections.Counter)
with open(OUT[0], "w") as o:
    o.write("stem\tcohort\tbarrnap_hit\tlength\tsilva_target\tidentity\tn_tied_top_hits\tgroup\ttaxonomy\n")
    for h in q:
        s = h.split("|", 1)[0]; c = coh[s]; f = top[h]; g, tax = group(f[1])
        idt = f[2] if f[1] != "*" else "NA"
        o.write("%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\n" % (s, c, h.split("|", 1)[1], q[h], f[1].split(" ")[0], idt, ties[h], g, tax))
        for k in (c, "all"): hits[k][g] += 1; samp[k][g].add(s)
        if idt != "NA": ids[g].append(float(idt))
        if g != "no hit": lowest[g][tax.split(";")[-1]] += 1
C = ["WF22", "WF23", "WF24", "all"]; NS = {c: sum(1 for s in coh if c == "all" or coh[s] == c) for c in C}
with open(OUT[1], "w") as o:
    o.write("group\t" + "\t".join("n_18S_%s\tpct_18S_%s\tn_samples_%s" % (c, c, c) for c in C) + "\n")
    for g in GROUPS + ["total"]:
        cells = []
        for c in C:
            n = sum(hits[c].values()) if g == "total" else hits[c][g]; tot = sum(hits[c].values())
            ns = len(set().union(*samp[c].values())) if g == "total" else len(samp[c][g])
            cells.append("%d\t%s\t%d" % (n, "%.2f" % (100.0 * n / tot) if tot else "NA", ns))
        o.write(g + "\t" + "\t".join(cells) + "\n")
print("18S hits >= 500 bp classified: %d (samples with any: %d of 88)" % (len(q), len(set().union(*samp["all"].values()))))
print("%-24s " % "group" + " ".join("%22s" % ("%s (n=%d samples)" % (c, NS[c])) for c in C))
print("%-24s " % "" + " ".join("%22s" % "18S hits (%) | samples" for c in C))
for g in GROUPS + ["total"]:
    row = []
    for c in C:
        n = sum(hits[c].values()) if g == "total" else hits[c][g]; tot = sum(hits[c].values())
        ns = len(set().union(*samp[c].values())) if g == "total" else len(samp[c][g])
        row.append("%22s" % ("%d (%.1f%%) | %d" % (n, 100.0 * n / tot if tot else 0, ns)))
    print("%-24s " % g + " ".join(row))
print("\ntop-hit identity by group (median, 10th pct) and the most common lowest SILVA taxa:")
for g in GROUPS[:-1]:
    v = sorted(ids[g])
    if not v: continue
    print("  %-24s n=%d median %.1f%% p10 %.1f%% | %s" % (g, len(v), v[len(v) // 2], v[len(v) // 10], "; ".join("%s %d" % x for x in lowest[g].most_common(5))))
for p in OUT: print("wrote", p)
