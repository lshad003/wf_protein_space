# Figure 1: a, supported representatives by annotation class and contig origin (alluvial); b, accumulation of strict
# microbial representatives over metagenomes (presence = a member on a strict-eligible contig; dashed: any member);
# c, strict representatives by cohort combination, stacked by class; d, per-metagenome share of catalog-mapped reads
# assigned to strict representatives. Inputs: results/step34a_origin_by_class.tsv, results/fig_prep_accum_v2.tsv,
# results/fig_prep_cohort_sharing.tsv, results/fig_prep_read_retention.tsv (Step 34k). Output: figures/fig1_catalog.*
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
fs.setup(); R = fs.W + "/results/"
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
CL = ["K", "KWP", "GU", "EU"]; CC = {"K": "#4477AA", "KWP": "#66CCEE", "GU": "#CCBB44", "EU": "#AA3377"}
OR = ["Tiara-prokaryotic", "<1 kb unclassified", "Tiara-unknown", "Tiara-eukaryotic other", "fly", "frog"]
OC = dict(zip(OR, ["#0072B2", "#BBBBBB", "#888888", "#E69F00", "#D55E00", "#CC79A7"]))
CO = {"pooled": "#222222", "WF22": "#332288", "WF23": "#88CCEE", "WF24": "#117733"}
oc = rows("step34a_origin_by_class.tsv"); use = [r for r in oc if r.get("len_bin") == "all"] or oc
F = {(c, o): sum(int(r["n_" + c]) for r in use if r["origin"] == o) for c in CL for o in OR}
LT = {c: sum(F[(c, o)] for o in OR) for c in CL}; RT = {o: sum(F[(c, o)] for c in CL) for o in OR}; TOT = sum(LT.values())
print("supported reps:", TOT)
fig = plt.figure(figsize=(7.08, 6.0))
g = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], width_ratios=[1.2, 0.8, 0.6], hspace=0.45, wspace=0.6)
ax = fig.add_subplot(g[0, 0]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
gap = 0.012; sc = (1 - gap * 5) / TOT; xL, xR, wN = 0.18, 0.52, 0.04
def stack(keys, tot):
    y, pos = 1.0, {}
    for k in keys: h = tot[k] * sc; pos[k] = [y - h, y]; y -= h + gap
    return pos
PL = stack(CL, LT); PR = stack(OR, RT); curL = {c: PL[c][1] for c in CL}; curR = {o: PR[o][1] for o in OR}
for c in CL:
    for o in OR:
        h = F[(c, o)] * sc
        if h <= 0: continue
        a0, a1, b0, b1 = curL[c] - h, curL[c], curR[o] - h, curR[o]; curL[c] -= h; curR[o] -= h
        x0, x1 = xL + wN, xR; xm = (x0 + x1) / 2
        verts = [(x0, a1), (xm, a1), (xm, b1), (x1, b1), (x1, b0), (xm, b0), (xm, a0), (x0, a0), (x0, a1)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=CC[c], alpha=0.45, linewidth=0))
for c in CL:
    ax.add_patch(Rectangle((xL, PL[c][0]), wN, PL[c][1] - PL[c][0], color=CC[c], linewidth=0))
    ax.text(xL - 0.015, sum(PL[c]) / 2, "%s\n%s" % (c, format(LT[c], ",")), ha="right", va="center", fontsize=5.5)
for o in OR:
    ax.add_patch(Rectangle((xR, PR[o][0]), wN, PR[o][1] - PR[o][0], color=OC[o], linewidth=0))
    ax.text(xR + wN + 0.015, sum(PR[o]) / 2, "%s  %s" % (o, format(RT[o], ",")), ha="left", va="center", fontsize=5.5)
ax.text(xL + wN / 2, 1.04, "annotation class", ha="center", fontsize=6); ax.text(xR + wN / 2, 1.04, "contig origin", ha="center", fontsize=6)
fs.letter(ax, "a")
ac = rows("fig_prep_accum_v2.tsv"); ax = fig.add_subplot(g[0, 1:])
for sn in ("WF23", "WF24", "WF22", "pooled"):
    d = [r for r in ac if r["rule"] == "strict_members" and r["set"] == sn]; x = [int(r["n_metagenomes"]) for r in d]
    m, lo, hi = (np.array([float(r[c]) for r in d]) / 1e6 for c in ("mean", "q025", "q975"))
    ax.fill_between(x, lo, hi, color=CO[sn], alpha=0.2, linewidth=0)
    ax.plot(x, m, color=CO[sn], linewidth=1.1 if sn == "pooled" else 0.8, label="%s (%d)" % ("all cohorts" if sn == "pooled" else sn, len(x)))
    ax.text(x[-1] + 1, m[-1], "%.2fM" % m[-1], fontsize=5, va="center", color=CO[sn])
d = [r for r in ac if r["rule"] == "all_members" and r["set"] == "pooled"]
ax.plot([int(r["n_metagenomes"]) for r in d], [float(r["mean"]) / 1e6 for r in d], color="#222222", linestyle=(0, (2, 2)), linewidth=0.7, label="all cohorts, any member")
ax.set_xlabel("metagenomes added (200 random orders)"); ax.set_ylabel("microbial representatives (millions)")
ax.set_xlim(0, 100); ax.set_ylim(bottom=0); ax.legend(frameon=False, fontsize=5, loc="lower right"); fs.letter(ax, "b")
sh = [r for r in rows("fig_prep_cohort_sharing.tsv") if r["rule"] == "strict_members"]
COMB = ["WF22", "WF23", "WF24", "WF22+WF23", "WF22+WF24", "WF23+WF24", "all3"]
LAB = ["WF22\nonly", "WF23\nonly", "WF24\nonly", "WF22+\nWF23", "WF22+\nWF24", "WF23+\nWF24", "all\nthree"]
ax = fig.add_subplot(g[1, 0:2]); bot = np.zeros(len(COMB)); eu = np.zeros(len(COMB))
for c in CL:
    v = np.array([sum(int(r["n_reps"]) for r in sh if r["combination"] == k and r["class"] == c) for k in COMB], dtype=float)
    ax.bar(range(len(COMB)), v / 1e3, bottom=bot / 1e3, color=CC[c], width=0.65, linewidth=0, label=c); bot += v
    if c == "EU": eu = v
for i in range(len(COMB)):
    ax.text(i, bot[i] / 1e3 + 6, "%s\nEU %.1f%%" % (format(int(bot[i]), ","), 100 * eu[i] / bot[i]), ha="center", va="bottom", fontsize=4.8)
print("strict reps with no eligible member:", sum(int(r["n_reps"]) for r in sh if r["combination"] == "none"))
ax.set_xticks(range(len(COMB))); ax.set_xticklabels(LAB, fontsize=5.5); ax.set_ylabel("microbial representatives (thousands)")
ax.set_ylim(0, bot.max() / 1e3 * 1.18); ax.legend(frameon=False, fontsize=5, loc="upper right", ncol=4); fs.letter(ax, "c")
rr = rows("fig_prep_read_retention.tsv"); ax = fig.add_subplot(g[1, 2]); rng = np.random.default_rng(1)
for i, c in enumerate(("WF22", "WF23", "WF24")):
    sel = [r for r in rr if r["cohort"] == c]; v = [float(r["pct_retained"]) for r in sel]
    ax.boxplot([v], positions=[i], widths=0.5, showfliers=False, medianprops=dict(color="black", linewidth=0.8),
               boxprops=dict(linewidth=0.6), whiskerprops=dict(linewidth=0.6), capprops=dict(linewidth=0.6))
    ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=5, color=CO[c], alpha=0.85, zorder=3, linewidths=0)
    pooled = 100 * sum(float(r["strict_reads"]) for r in sel) / sum(float(r["catalog_reads"]) for r in sel)
    ax.text(i, 99.2, "%.1f%%" % pooled, ha="center", fontsize=5); print("%s pooled retention %.2f%%" % (c, pooled))
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["WF22", "WF23", "WF24"]); ax.set_ylim(80, 100)
ax.set_ylabel("catalog-mapped reads in microbial catalog (%)"); fs.letter(ax, "d")
fs.save(fig, "fig1_catalog")
