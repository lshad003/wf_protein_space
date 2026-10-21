# Figure 5: host factors. a-b, PCoA of WF22 metagenomes (Bray-Curtis on Hellinger-transformed abundances of genes in
# >= 35 of 44 samples, as Step 34e) for all microbial genes and for EU genes; colour = egg mass, marker = month,
# outline = sequencing run, lines join each froglet's months. c, PERMANOVA R2 for all genes vs EU genes (filled:
# p < 0.05). Inputs: results/fig_pcoa_wf22.tsv, results/batch_table_88.tsv, results/step34e_summary.tsv.
# Output: figures/fig5_host_factors.pdf/.png
# Optional arguments (Step 34s): <dir> <suffix>: read-based inputs (fig_pcoa_wf22.tsv, step34e_summary.tsv) are read from <dir>; output
#   figures/fig5_host_factors<suffix>.pdf/.png. Default (no arguments) unchanged.
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
fs.setup(); R = fs.W + "/results/"
def rows(f, d=None): return list(csv.DictReader(open((d or R) + f), delimiter="\t"))
D = sys.argv[1].rstrip("/") + "/" if len(sys.argv) > 1 else R; TAG = sys.argv[2] if len(sys.argv) > 2 else ""
print("read-based inputs from", D)
def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None
pc = rows("fig_pcoa_wf22.tsv", D); bt = rows("batch_table_88.tsv"); bh = list(bt[0].keys())
sc = next(c for c in bh if c.lower() in ("stem", "sample", "sample_id"))
rc = next((c for c in bh if "run" in c.lower()), None) or next(c for c in bh if "batch" in c.lower())
print("batch table columns:", bh, "| using", sc, "and", rc); run = {r[sc]: r[rc] for r in bt}
EMC = {"1": "#1B9E77", "2": "#D95F02", "3": "#7570B3"}; MK = {"1": "o", "2": "s", "3": "^"}
runs = sorted({run[r["stem"]] for r in pc})
RC = dict(zip(runs, ["#000000", "#E69F00", "#56B4E9", "#CC79A7", "#999999", "#009E73"]))
print("WF22 metagenomes per run:", {k: sum(1 for r in pc if r["set"] == "strict_final" and run[r["stem"]] == k) for k in runs})
fig = plt.figure(figsize=(7.08, 5.0)); g = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.5, wspace=0.4)
for k, (lab, title) in enumerate((("strict_final", "all microbial genes"), ("strict_final_EU", "EU genes"))):
    ax = fig.add_subplot(g[0, k]); d = [r for r in pc if r["set"] == lab]
    for an in sorted({r["animal"] for r in d}):
        s = sorted([r for r in d if r["animal"] == an], key=lambda r: int(r["month"]))
        ax.plot([float(r["PC1"]) for r in s], [float(r["PC2"]) for r in s], color=EMC[str(s[0]["egg_mass"])], alpha=0.35, linewidth=0.5, zorder=1)
    for r in d:
        ax.scatter(float(r["PC1"]), float(r["PC2"]), s=20, marker=MK[str(r["month"])], color=EMC[str(r["egg_mass"])],
                   edgecolors=RC[run[r["stem"]]], linewidths=0.9, zorder=2)
    ax.set_xlabel("PCo1 (%.1f%%)" % float(d[0]["pct1"])); ax.set_ylabel("PCo2 (%.1f%%)" % float(d[0]["pct2"]))
    ax.text(0.02, 0.98, title, transform=ax.transAxes, fontsize=5.5, va="top"); fs.letter(ax, "ab"[k])
h = [Line2D([], [], marker="o", linestyle="", color=EMC[e], markersize=4, label="egg mass %s" % e) for e in "123"]
h += [Line2D([], [], marker=MK[m], linestyle="", color="#BBBBBB", markersize=4, label="month %s" % m) for m in "123"]
h += [Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=RC[x], markeredgewidth=1.0, markersize=4, label=x) for x in runs]
fig.legend(handles=h, ncol=6, frameon=False, fontsize=4.8, loc="upper center", bbox_to_anchor=(0.5, 1.03))
sm = rows("step34e_summary.tsv", D)
pr = [r for r in sm if r["kind"] == "permanova" and (r["term_or_coef"] in ("egg_mass", "month", "treatment") or (r["term_or_coef"] == "Model" and r["cohort"] in ("WF23", "WF24_n34")))
      and not any(k in r["test"].lower() for k in ("disp", "krus", "rich", "35"))
      and num(r["strict_final"]) is not None and num(r["strict_final_EU"]) is not None]
def lbl(r):
    t = r["test"].lower(); ts = "whole animal" if "whole" in t else "centroids" if "centroid" in t else "within animal" if "within" in t else "one sample per animal"
    m = "detection" if ("jacc" in r["metric_or_model"].lower() or "presence" in r["metric_or_model"].lower()) else "abundance"
    return "%s %s: %s (%s)" % (r["cohort"].replace("_n34", ""), m, "treatment" if r["term_or_coef"] == "Model" else r["term_or_coef"].replace("_", " "), ts)
ax = fig.add_subplot(g[1, :]); y = np.arange(len(pr))[::-1]
for yy, r in zip(y, pr):
    a, e, pa, pe = num(r["strict_final"]), num(r["strict_final_EU"]), num(r["strict_final_p"]), num(r["strict_final_EU_p"])
    ax.plot([a, e], [yy, yy], color="#BBBBBB", linewidth=0.8, zorder=1)
    ax.scatter(a, yy, s=16, color="#4477AA" if pa is not None and pa < 0.05 else "white", edgecolors="#4477AA", linewidths=0.8, zorder=2)
    ax.scatter(e, yy, s=16, color="#AA3377" if pe is not None and pe < 0.05 else "white", edgecolors="#AA3377", linewidths=0.8, zorder=2)
    print("R2 %-55s all %.3f (p %s) | EU %.3f (p %s)" % (lbl(r), a, r["strict_final_p"], e, r["strict_final_EU_p"]))
ax.set_yticks(y); ax.set_yticklabels([lbl(r) for r in pr], fontsize=5); ax.set_xlabel("PERMANOVA R\u00b2"); ax.set_xlim(left=0)
ax.legend(handles=[Line2D([], [], marker="o", linestyle="", color="#4477AA", markersize=4, label="all microbial genes"),
                   Line2D([], [], marker="o", linestyle="", color="#AA3377", markersize=4, label="EU genes"),
                   Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="#555555", markersize=4, label="p \u2265 0.05")],
          frameon=False, fontsize=5, loc="lower right"); fs.letter(ax, "c")
fs.save(fig, "fig5_host_factors" + TAG)
