# Step 33e: Tiara-eukaryotic contigs against read-level Kraken2 and contig-level UniRef50 LCA. Descriptive only:
# cohort is confounded with sequencing run.
# 1. Per sample (88, metadata/samples_88_treatment.tsv): Tiara-eukaryotic fraction of Tiara-classified contigs >= 3 kb
#    (eukarya + organelle over all Tiara contigs >= 3 kb; also by bp), next to Kraken2 PlusPFP read percentages from
#    results_kraken2/<stem>/<stem>.pluspfp.kraken2.report: unclassified (U), Eukaryota (D 2759), Chordata (P 7711),
#    Homo sapiens (S 9606), Streptophyta (P 35493). Percent = clade reads (col 2) / (unclassified + root reads).
#    Spearman (scipy) overall and within each cohort. Output results/step33e_tiara_vs_kraken.tsv.
# 2. Contigs >= 3 kb with a Tiara call and a domain in results/gene_taxonomy_n88.tsv (step20c: the domain of the
#    contig's UniRef50 LCA, reached through a supported rep on it; "unclassified" rows have no domain). One row per
#    contig; aborts if genes on one contig disagree. Cross-tab Tiara group (prok = bacteria + archaea + prokarya,
#    euk = eukarya + organelle, unknown) x domain (Bacteria / Archaea / Eukaryota / Viruses). Agreement = prok with
#    Bacteria or Archaea, euk with Eukaryota, over contigs with a prok/euk call and a B/A/E domain.
#    Output results/step33e_tiara_vs_lca.tsv.
import os, re, collections
from scipy.stats import spearmanr
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
K = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_kraken2"
OUT = [R + "step33e_tiara_vs_kraken.tsv", R + "step33e_tiara_vs_lca.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
GRP = {"bacteria": "prok", "archaea": "prok", "prokarya": "prok", "eukarya": "euk", "organelle": "euk", "unknown": "unknown"}
TAX = (("unclassified", "U", "0"), ("Eukaryota", "D", "2759"), ("Chordata", "P", "7711"), ("Homo_sapiens", "S", "9606"), ("Streptophyta", "P", "35493"))
coh = {}
with open(W + "/metadata/samples_88_treatment.tsv") as fh:
    next(fh)
    for l in fh: f = l.split("\t"); coh[f[0]] = f[1]
if len(coh) != 88: raise SystemExit("ABORT: %d stems" % len(coh))
tia, rows = {}, []
for s in coh:
    n = collections.Counter(); bp = collections.Counter()
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 2); ln = int(re.search(r"len=(\d+)", f[0]).group(1)); g = GRP[f[1]]
            if ln >= 3000: tia[(s, f[0].split()[0])] = g; n[g] += 1; bp[g] += ln
    rep = "%s/%s/%s.pluspfp.kraken2.report" % (K, s, s); got = {}
    for l in open(rep):
        f = l.rstrip("\n").split("\t")
        for nm, rk, tx in TAX + (("root", "R", "1"),):
            if f[3] == rk and f[4] == tx: got[nm] = int(f[1])
    miss = [nm for nm, _, _ in TAX + (("root", "R", "1"),) if nm not in got and nm not in ("Homo_sapiens", "Streptophyta", "Chordata")]
    if miss: raise SystemExit("ABORT: %s lacks %s" % (rep, miss))
    tot = got["unclassified"] + got["root"]
    rows.append([s, coh[s], sum(n.values()), n["euk"], n["euk"] / sum(n.values()), bp["euk"] / sum(bp.values())]
                + [100.0 * got.get(nm, 0) / tot for nm, _, _ in TAX] + [tot])
with open(OUT[0], "w") as o:
    o.write("stem\tcohort\tn_tiara_contigs_ge3kb\tn_euk_ge3kb\tfrac_euk_contigs_ge3kb\tfrac_euk_bp_ge3kb\t"
            + "\t".join("kraken_pct_" + nm for nm, _, _ in TAX) + "\tkraken_reads_total\n")
    for r in rows: o.write("\t".join([r[0], r[1], str(r[2]), str(r[3]), "%.4f" % r[4], "%.4f" % r[5]] + ["%.4f" % x for x in r[6:11]] + [str(r[11])]) + "\n")
def med(v): v = sorted(v); return v[len(v) // 2]
print("1. per sample (n=88): Tiara-eukaryotic fraction of contigs >= 3 kb vs Kraken2 PlusPFP read %")
for c in ("all", "WF22", "WF23", "WF24"):
    s = [r for r in rows if c == "all" or r[1] == c]
    print("  %s n=%d | Tiara euk frac median %.3f (range %.3f-%.3f) | Kraken unclassified%% median %.2f (range %.2f-%.2f) | Eukaryota%% median %.2f"
          % (c, len(s), med([r[4] for r in s]), min(r[4] for r in s), max(r[4] for r in s), med([r[6] for r in s]), min(r[6] for r in s), max(r[6] for r in s), med([r[7] for r in s])))
print("  Spearman rho (p) of Tiara euk fraction (contigs; bp in brackets) against:")
for i, (nm, _, _) in enumerate(TAX):
    out = []
    for c in ("all", "WF22", "WF23", "WF24"):
        s = [r for r in rows if c == "all" or r[1] == c]
        a = spearmanr([r[4] for r in s], [r[6 + i] for r in s]); b = spearmanr([r[5] for r in s], [r[6 + i] for r in s])
        out.append("%s n=%d %.2f (%.2g) [%.2f]" % (c, len(s), a[0], a[1], b[0]))
    print("    %-13s " % nm + " | ".join(out))
dom, bad = {}, 0
with open(R + "gene_taxonomy_n88.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); d = f[4]
        if d == "unclassified": continue
        s, rest = f[0].split("__", 1); k = (s, rest[:rest.rfind("_")])
        if k in dom and dom[k] != d: bad += 1
        dom[k] = d
if bad: raise SystemExit("ABORT: %d genes disagree with another gene's domain on the same contig" % bad)
xt = collections.Counter(); both = [k for k in dom if k in tia]
for k in both: xt[(tia[k], dom[k])] += 1
G, D = ("prok", "euk", "unknown"), ("Bacteria", "Archaea", "Eukaryota", "Viruses")
with open(OUT[1], "w") as o:
    o.write("tiara_group\t" + "\t".join(D) + "\ttotal\n")
    for g in G: o.write(g + "\t" + "\t".join(str(xt[(g, d)]) for d in D) + "\t%d\n" % sum(xt[(g, d)] for d in D))
    o.write("total\t" + "\t".join(str(sum(xt[(g, d)] for g in G)) for d in D) + "\t%d\n" % len(both))
print("\n2. contigs >= 3 kb with a Tiara call and a UniRef50-LCA domain (gene_taxonomy_n88): %d (of %d contigs with a domain, any length)" % (len(both), len(dom)))
print("   %-8s %10s %10s %10s %10s" % (("",) + D))
for g in G: print("   %-8s %10d %10d %10d %10d" % ((g,) + tuple(xt[(g, d)] for d in D)))
ag = xt[("prok", "Bacteria")] + xt[("prok", "Archaea")] + xt[("euk", "Eukaryota")]
den = sum(xt[(g, d)] for g in ("prok", "euk") for d in ("Bacteria", "Archaea", "Eukaryota"))
pe = sum(xt[("prok", d)] for d in ("Bacteria", "Archaea")); ee = sum(xt[("euk", d)] for d in ("Bacteria", "Archaea", "Eukaryota"))
print("   agreement (prok-Bacteria/Archaea, euk-Eukaryota) over prok/euk x B/A/E: %d/%d (%.2f%%)" % (ag, den, 100.0 * ag / den))
print("   Tiara-prok with a B/A/E domain that is Bacteria/Archaea: %.2f%% | Tiara-euk with a B/A/E domain that is Eukaryota: %d/%d (%.2f%%)"
      % (100.0 * pe / sum(xt[("prok", d)] for d in ("Bacteria", "Archaea", "Eukaryota")), xt[("euk", "Eukaryota")], ee, 100.0 * xt[("euk", "Eukaryota")] / ee))
print("   LCA-Eukaryota contigs called euk by Tiara: %d/%d (%.2f%%)" % (xt[("euk", "Eukaryota")], sum(xt[(g, "Eukaryota")] for g in G), 100.0 * xt[("euk", "Eukaryota")] / sum(xt[(g, "Eukaryota")] for g in G)))
for p in OUT: print("wrote", p)
