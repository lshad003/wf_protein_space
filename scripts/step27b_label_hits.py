# Step 27b: label Basidiobolus-vs-catalog hits (step27a) by eggNOG lineage, contig domain, source stem and treatment.
# Hits: results/basid_vs_catalog_hits.tsv (query = strain|protein, target = catalog 95% rep, 89-sample rep fasta;
#   mmseqs search -e 1e-5 --min-seq-id 0.9 -c 0.8 --cov-mode 1, i.e. >= 90% identity, >= 80% target coverage).
#   Columns: query,target,pident,alnlen,qlen,tlen,qcov,tcov,evalue,bits.
# eggNOG (Step 8): results/eggnog2/LsPS.{0..619}.emapper.annotations; eggNOG_OGs (col 5), max_annot_lvl (col 6).
#   Lineage from the taxonomic levels in eggNOG_OGs: 2|Bacteria -> Bacteria; 4751|Fungi -> Fungi;
#   2759|Eukaryota -> other Eukaryota; 2157|Archaea -> Archaea; annotated with root only -> root_only;
#   no row -> none. Bacteria together with a eukaryotic level is counted as a conflict and labelled "conflict".
# Contig domain: results/gene_taxonomy_n88.tsv (rep, domain), covers 48 stems only; NA where absent.
# Stem = target before "__". Cohort and treatment: metadata/samples_88_treatment.tsv; WF24 codes resolved with
#   metadata/wf24_treatment_key.tsv (code 4 = Control). Stems not in the 88 (UHM586.41010, excluded contaminated
#   WF24 control, metadata/wf24_excluded.tsv) get treatment "excluded" and count as neither inoculated nor control.
# Supported: results/supported_reps_95_n88.txt. In n88 95% tier: results/reps_95_n88.txt.
# Bacterial target: eggNOG lineage Bacteria OR contig domain Bacteria.
# Best hit per query: highest bits, then highest pident, then first in file.
# Read-only on every _n88 file.
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results"
M = W + "/metadata"
OUT1 = R + "/basid_vs_catalog_labelled.tsv"
OUT2 = R + "/basid_bacterial_candidates.tsv"
for o in (OUT1, OUT2):
    if os.path.exists(o): raise SystemExit("REFUSING: " + o + " exists")

hits = [l.rstrip("\n").split("\t") for l in open(R + "/basid_vs_catalog_hits.tsv")]
tg = {h[1] for h in hits}
print("hit rows:", len(hits), "| queries:", len({h[0] for h in hits}), "| targets:", len(tg))

# eggNOG lineage
LV = {"2": "Bacteria", "2157": "Archaea", "2759": "Eukaryota", "4751": "Fungi"}
egg = {}
nfiles = 0
for n in range(620):
    p = R + "/eggnog2/LsPS.%d.emapper.annotations" % n
    if not os.path.exists(p): continue
    nfiles += 1
    for l in open(p):
        if l[0] == "#": continue
        k, rest = l.split("\t", 1)
        if k not in tg: continue
        f = l.rstrip("\n").split("\t")
        lv = {x.split("@", 1)[1].split("|", 1)[0] for x in f[4].split(",") if "@" in x}
        b, e = "2" in lv, ("2759" in lv or "4751" in lv)
        if b and e: lin = "conflict"
        elif b: lin = "Bacteria"
        elif "4751" in lv: lin = "Fungi"
        elif "2759" in lv: lin = "other Eukaryota"
        elif "2157" in lv: lin = "Archaea"
        else: lin = "root_only"
        egg[k] = (lin, f[5])
print("eggNOG annotation files read:", nfiles)
if nfiles != 619: raise SystemExit("ABORT: expected 619 eggNOG files")

dom = {}
for i, l in enumerate(open(R + "/gene_taxonomy_n88.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    if f[0] in tg: dom[f[0]] = f[4]
sup = {l.strip() for l in open(R + "/supported_reps_95_n88.txt") if l.strip() in tg}
rep88 = {l.strip() for l in open(R + "/reps_95_n88.txt") if l.strip() in tg}

key = {}
for i, l in enumerate(open(M + "/wf24_treatment_key.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t"); key[f[0]] = f[1]
meta = {}
for i, l in enumerate(open(M + "/samples_88_treatment.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    t = key[f[2]] if f[1] == "WF24" else f[2]
    meta[f[0]] = (f[1], t, "control" if t.lower() == "control" else "inoculated")
for i, l in enumerate(open(M + "/wf24_excluded.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    meta.setdefault(f[0], ("WF24", "excluded", "excluded"))

def lab(t):
    stem = t.split("__", 1)[0]
    if stem not in meta: raise SystemExit("ABORT: stem without metadata: " + stem)
    lin, mal = egg.get(t, ("none", "NA"))
    d = dom.get(t, "NA")
    bac = "yes" if (lin == "Bacteria" or d == "Bacteria") else "no"
    c, tr, grp = meta[stem]
    return [lin, mal, d, bac, stem, c, tr, grp, "yes" if t in sup else "no", "yes" if t in rep88 else "no"]

best = {}
for i, h in enumerate(hits):
    s = (float(h[9]), float(h[2]), -i)
    if h[0] not in best or s > best[h[0]][0]: best[h[0]] = (s, i)
bi = {v[1] for v in best.values()}

HDR = ["query", "strain", "target", "pident", "alnlen", "qlen", "tlen", "qcov", "tcov", "evalue", "bits",
       "egg_lineage", "max_annot_lvl", "contig_domain", "bacterial", "stem", "cohort", "treatment", "group",
       "supported", "in_n88_95_tier", "is_best_hit"]
rows = []
for i, h in enumerate(hits):
    rows.append([h[0], h[0].split("|", 1)[0]] + h[1:] + lab(h[1]) + ["yes" if i in bi else "no"])
with open(OUT1, "w") as fh:
    fh.write("\t".join(HDR) + "\n")
    for r in rows: fh.write("\t".join(r) + "\n")

C = {c: i for i, c in enumerate(HDR)}
bb = [rows[i] for i in sorted(bi) if rows[i][C["bacterial"]] == "yes"]
bb.sort(key=lambda r: (-float(r[C["pident"]]), -float(r[C["tcov"]]), -float(r[C["bits"]])))
with open(OUT2, "w") as fh:
    fh.write("\t".join(HDR[:-1]) + "\n")
    for r in bb: fh.write("\t".join(r[:-1]) + "\n")

# summary
tl = collections.Counter(lab(t)[0] for t in tg)
td = collections.Counter(lab(t)[2] for t in tg)
print("\ntargets by eggNOG lineage:", dict(tl))
print("targets by contig domain:", dict(td))
print("targets supported:", len(sup), "of", len(tg), "| in n88 95% tier:", len(rep88))
print("targets bacterial (egg OR contig):", sum(lab(t)[3] == "yes" for t in tg),
      "| egg only:", sum(lab(t)[0] == "Bacteria" and lab(t)[2] != "Bacteria" for t in tg),
      "| contig only:", sum(lab(t)[0] != "Bacteria" and lab(t)[2] == "Bacteria" for t in tg),
      "| both:", sum(lab(t)[0] == "Bacteria" and lab(t)[2] == "Bacteria" for t in tg))
nq = len(best)
print("\nqueries with any hit:", nq)
for s, n in sorted(collections.Counter(q.split("|", 1)[0] for q in best).items()): print("  ", s, n)
print("queries whose best hit is bacterial:", len(bb), "(%.2f%% of queries with a hit)" % (100.0 * len(bb) / nq))
print("  of which best-hit target supported:", sum(r[C["supported"]] == "yes" for r in bb))
for s, n in sorted(collections.Counter(r[C["strain"]] for r in bb).items()): print("  ", s, n)
print("best-hit bacterial by source frog group (stem of the target representative):")
for g in ("inoculated", "control", "excluded"):
    a = [r for r in bb if r[C["group"]] == g]
    print("  %-10s all %d | supported %d" % (g, len(a), sum(r[C["supported"]] == "yes" for r in a)))
print("best-hit bacterial by cohort and treatment:")
for (c, t), n in sorted(collections.Counter((r[C["cohort"]], r[C["treatment"]]) for r in bb).items()):
    print("  ", c, t, n)
print("\ntop 30 bacterial best hits (pident desc, tcov desc, bits desc):")
S = ["query", "target", "pident", "tcov", "egg_lineage", "contig_domain", "stem", "treatment"]
print("\t".join(S))
for r in bb[:30]: print("\t".join(r[C[c]] for c in S))
print("\nwrote", OUT1, "and", OUT2)
