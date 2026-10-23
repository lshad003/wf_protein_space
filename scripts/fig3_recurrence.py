# Figure 3: recurrence at three levels (microbial catalog, EU), six panels at 180 mm width, one font size (6 pt).
# a, level 1: 30% families by animals and cohorts, recurrent zone outlined; b, level 2: within-year recurrent EU and K reps
# across years (UpSet); c, level 3: WF22 egg-mass-restricted EU reps across the 15 WF22 froglets; d, egg-mass-restricted
# EU and K reps per WF22 egg mass vs the 95th percentile of each class's null; e, WF22 reps meeting the egg-mass rule in
# 1, 2 or 3 egg masses; f, tank-restricted EU and K reps (WF24 EM4, WF23 EM3). Panels d-f from scripts/fig3_tank_eggmass.py.
# Inputs: results/step34c2c_strict_v2_families.tsv, fig_prep_year_flags.tsv, fig_prep_wf22_em_presence.tsv,
# results/step34t/ (wf22_eggmass_nulls.tsv, wf22_eggmass_breadth.tsv, tank_run_check.tsv, tank_run_crosstab.tsv).
# Output: figures/fig3_recurrence.*
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import fig3_tank_eggmass as P
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import Rectangle, Patch
fs.setup(); R = fs.W + "/results/"
FS = 6
plt.rcParams.update({"font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
                     "legend.fontsize": FS, "legend.title_fontsize": FS})
def letter(ax, s, x=-0.12):
    ax.text(x, 1.06, s, transform=ax.transAxes, fontsize=FS, fontweight="bold", va="bottom", ha="left")
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
EU, KC = P.EU, P.KC; EMC = {"1": "#1B9E77", "2": "#D95F02", "3": "#7570B3"}
TRC = {"control": "#E8E8E8", "stp1710.7": "#333333", "stp1717.1": "#999999"}
fig = plt.figure(figsize=(180 / 25.4, 8.6))
G = fig.add_gridspec(4, 1, height_ratios=[1, 1.15, 0.16, 1], hspace=0.42)   # row 3 is a spacer for the panel c legend
r1 = G[0].subgridspec(1, 2, width_ratios=[1, 1.25], wspace=0.38)
r2 = G[1].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.32)
r3 = G[3].subgridspec(1, 2, width_ratios=[0.75, 2.0], wspace=0.55)
# a
fam = [r for r in rows("step34c2c_strict_v2_families.tsv") if r["excluded"] == "0"]
an = np.array([int(r["n_animals"]) for r in fam]); co = np.array([int(r["n_cohorts"]) for r in fam])
AB = [(1, 1, "1"), (2, 4, "2-4"), (5, 9, "5-9"), (10, 19, "10-19"), (20, 39, "20-39"), (40, 999, "≥40")]
H = np.array([[((an >= lo) & (an <= hi) & (co == k)).sum() for lo, hi, _ in AB] for k in (1, 2, 3)])
nrec = sum(1 for r in fam if r["recurrent"] == "1"); print("families %d, recurrent %d, outlined cells %d" % (len(fam), nrec, int(H[1:, 3:].sum())))
ax = fig.add_subplot(r1[0])
ax.imshow(np.where(H > 0, H, np.nan), cmap="RdPu", norm=mc.LogNorm(vmin=1, vmax=H.max()), aspect="auto", origin="lower")
for i in range(3):
    for j in range(len(AB)):
        if H[i, j]: ax.text(j, i, format(int(H[i, j]), ","), ha="center", va="center", fontsize=FS, color="white" if H[i, j] > 3000 else "black")
ax.add_patch(Rectangle((2.5, 0.5), 3, 2, fill=False, edgecolor="black", linewidth=1.0))
ax.set_xticks(range(len(AB))); ax.set_xticklabels([a[2] for a in AB]); ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["1", "2", "3"])
ax.set_xlabel("animals per family"); ax.set_ylabel("cohorts per family")
ax.text(4, 2.62, "recurrent (outlined): %s families" % format(nrec, ","), ha="center", fontsize=FS); letter(ax, "a")
# b
fl = rows("fig_prep_year_flags.tsv"); Y = ["WF22", "WF23", "WF24"]
COMB = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
cnt = {s: [sum(1 for r in fl if r["set"] == s and tuple(int(r["wr_" + y]) for y in Y) == c) for c in COMB] for s in ("EU", "K")}
for s in cnt: print("level 2 %s intersections:" % s, dict(zip(["22", "23", "24", "22+23", "22+24", "23+24", "all"], cnt[s])))
sg = r1[1].subgridspec(2, 1, height_ratios=[2.2, 1], hspace=0.05)
ax = fig.add_subplot(sg[0]); x = np.arange(len(COMB)); w = 0.38
ax.bar(x - w / 2, cnt["EU"], w, color=EU, linewidth=0, label="EU"); ax.bar(x + w / 2, cnt["K"], w, color=KC, linewidth=0, label="K (equal-size random set)")
for i in x:
    for dx, s in ((-w / 2, "EU"), (w / 2, "K")):
        ax.text(i + dx, cnt[s][i] + 40, format(cnt[s][i], ","), ha="center", va="bottom", fontsize=FS, rotation=90)
ax.set_xticks([]); ax.set_xlim(-0.6, len(COMB) - 0.4); ax.set_ylim(0, 1.22 * max(cnt["EU"] + cnt["K"]))
ax.set_ylabel("within-year recurrent\nrepresentatives"); ax.legend(frameon=False, loc="lower left", ncol=2, bbox_to_anchor=(0.0, 1.0)); letter(ax, "b")
ax2 = fig.add_subplot(sg[1], sharex=ax)
for i, c in enumerate(COMB):
    ys = [2 - j for j in range(3) if c[j]]
    ax2.scatter([i] * 3, [2, 1, 0], s=10, color="#DDDDDD", zorder=1); ax2.scatter([i] * len(ys), ys, s=10, color="#333333", zorder=2)
    if len(ys) > 1: ax2.plot([i, i], [min(ys), max(ys)], color="#333333", linewidth=0.8, zorder=2)
ax2.set_yticks([2, 1, 0]); ax2.set_yticklabels(Y); ax2.set_ylim(-0.5, 2.5); ax2.tick_params(axis="x", length=0)
for sp in ("bottom", "left"): ax2.spines[sp].set_visible(False)
# c
pe = list(csv.reader(open(R + "fig_prep_wf22_em_presence.tsv"), delimiter="\t")); hd = [h.split("|") for h in pe[0][2:]]; body = pe[1:]
cols = sorted(range(len(hd)), key=lambda i: (hd[i][1], hd[i][2].lower(), hd[i][0]))
body.sort(key=lambda r: (r[1], -sum(int(v) for v in r[2:]))); M = np.array([[int(r[2 + i]) for i in cols] for r in body])
sg = r2[0].subgridspec(2, 1, height_ratios=[0.09, 1], hspace=0.04); axs = fig.add_subplot(sg[0]); ax = fig.add_subplot(sg[1])
strip = np.array([[mc.to_rgb(EMC[hd[i][1]]) for i in cols], [mc.to_rgb(TRC.get(hd[i][2].lower(), "#FFFFFF")) for i in cols]])
axs.imshow(strip, aspect="auto"); axs.set_xticks([]); axs.set_yticks([0, 1]); axs.set_yticklabels(["egg mass", "treatment"])
for sp in axs.spines.values(): sp.set_visible(False)
letter(axs, "c", x=-0.2)
ax.imshow(M, aspect="auto", cmap=mc.ListedColormap(["#F2F2F2", EU]), interpolation="nearest", rasterized=True)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([hd[i][0] for i in cols], rotation=90); ax.set_yticks([])
start = 0
for e in ("EM1", "EM2", "EM3"):
    k = sum(1 for r in body if r[1] == e)
    if k: ax.axhline(start + k - 0.5, color="black", linewidth=0.4); ax.text(-0.7, start + k / 2, "%s %d" % (e, k), ha="right", va="center", fontsize=FS)
    start += k
for j in range(1, len(cols)):
    if hd[cols[j]][1] != hd[cols[j - 1]][1]: ax.axvline(j - 0.5, color="black", linewidth=0.4)
print("level 3 heatmap: %d reps x %d froglets" % M.shape)
hl = [Patch(color=EMC[e], label="EM%s" % e) for e in "123"] + [Patch(color=TRC[t], label=l) for t, l in (("control", "control"), ("stp1717.1", "STP1717.1"), ("stp1710.7", "STP1710.7"))]
ax.legend(handles=hl, ncol=3, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), handlelength=1, columnspacing=0.8)
# d
ax = fig.add_subplot(r2[1]); hl = P.eggmass_panel(ax, FS); letter(ax, "d")
ax.legend(handles=hl, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.92))
# e
ax = fig.add_subplot(r3[0]); P.breadth_panel(ax, FS); letter(ax, "e", x=-0.45)
# f
sf = r3[1].subgridspec(1, 2, wspace=0.12); a1 = fig.add_subplot(sf[0]); a2 = fig.add_subplot(sf[1], sharey=a1)
P.tank_panel([a1, a2], P.CON, P.MODES, FS); a2.tick_params(axis="y", labelleft=False); letter(a1, "f", x=-0.25)
fs.save(fig, "fig3_recurrence")
P.legend_d()
