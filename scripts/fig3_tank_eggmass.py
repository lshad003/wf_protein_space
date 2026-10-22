# Figure 3 additions (Step 34t): tank panel styled like Fig. 3d, observed tank-restricted EU and K reps vs permutation null
# (mean +/- sd; each class's own one-sided p under its point), WF24 egg mass 4 and WF23 egg mass 3, assembly presence and
# >= 5 reads at c_id95 (panel f); and a stacked bar of WF22 months 2-3 reps meeting the Step 34c3 rule in exactly 1, 2 or
# 3 egg masses (assembly; panel e).
# Inputs: results/step34t/tank_restriction_c_id95.tsv, results/step34t/wf22_eggmass_breadth.tsv.
# Outputs: figures/fig3_tank_panel.*, figures/fig3_wf22_eggmass_breadth.*
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
fs.setup(); R = fs.W + "/results/step34t/"
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
EU, KC = "#AA3377", "#4477AA"
T = rows("tank_restriction_c_id95.tsv")
CON = [("WF24_EM4_tank", "WF24 EM4"), ("WF23_EM3_tank", "WF23 EM3")]
MODES = [("assembly", "assembly presence"), ("reads>=5_c_id95", "≥ 5 reads (≥ 95% identity)")]
fig, axes = plt.subplots(1, 2, figsize=(5.2, 1.9), gridspec_kw=dict(wspace=0.75))
for ax, (mode, mlab) in zip(axes, MODES):
    xmax = 0
    for k, (s, col) in enumerate((("EU", EU), ("K", KC))):
        rr = [next(r for r in T if r["contrast"] == c and r["presence"] == mode and r["set"] == s) for c, _ in CON]
        y = np.arange(len(rr)) + (0.15 if k == 0 else -0.15)
        ax.errorbar([float(r["null_mean"]) for r in rr], y, xerr=[float(r["null_sd"]) for r in rr], fmt="|", color="#AAAAAA", elinewidth=1.5, capsize=0)
        ax.scatter([int(r["observed_restricted"]) for r in rr], y, s=18, color=col, zorder=3, label=s + " observed")
        xmax = max(xmax, max(int(r["observed_restricted"]) for r in rr))
        for yy, r in zip(y, rr):
            ax.text(int(r["observed_restricted"]) + 0.02 * xmax + 15, yy, format(int(r["observed_restricted"]), ","), va="center", fontsize=5)
            ax.text(int(r["observed_restricted"]), yy + 0.09, "p = %.3f" % float(r["p_one_sided"]), ha="center", va="top", fontsize=4.2, color=col)
    for i, (c, _) in enumerate(CON):
        d = next(r for r in T if r["contrast"] == c and r["presence"] == mode and r["set"] == "D_EU_minus_K")
        ax.text(1.0, i + 0.42, "EU − K = %s, p = %.3f" % (format(int(d["observed_restricted"]), ","), float(d["p_one_sided"])),
                transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=4.8)
    tl = []
    for c, lab in CON:
        r = next(r for r in T if r["contrast"] == c and r["presence"] == mode and r["set"] == "EU")
        tl.append("%s\n(%s frogs)" % (lab, r["animals_per_tank"].replace(";", ", ").replace(":", ": ")))
    ax.set_xlim(0, xmax * 1.25); ax.set_ylim(-0.5, 1.7); ax.set_yticks(range(len(CON))); ax.set_yticklabels(tl, fontsize=5); ax.invert_yaxis()
    ax.set_title(mlab, fontsize=6); ax.set_xlabel("tank-restricted representatives\n(grey: permutation null, mean ± sd)")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=5, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.12))
fig.text(0.5, -0.27, "p under each point: that class vs its own null; EU − K p: one-sided permutation test of the EU − K difference\n"
         "(both one-sided, 1,000 tank-label permutations within the egg mass)", ha="center", fontsize=5)
fs.letter(axes[0], "f")
fs.save(fig, "fig3_tank_panel")
B = rows("wf22_eggmass_breadth.tsv")
fig, ax = plt.subplots(figsize=(1.9, 2.0))
SH = {"1": "#DDDDDD", "2": "#999999", "3": "#333333"}
for i, (s, col) in enumerate((("EU", EU), ("K", KC))):
    bot = 0
    for k in ("1", "2", "3"):
        r = next(r for r in B if r["set"] == s and r["n_egg_masses_meeting_rule"] == k); n = int(r["n_reps"])
        ax.bar(i, n, 0.6, bottom=bot, color=SH[k], edgecolor=col, linewidth=1.0)
        if n: ax.text(i + 0.33, bot + n / 2, format(n, ","), va="center", fontsize=5)
        bot += n
ax.set_xticks([0, 1]); ax.set_xticklabels(["EU", "K"]); ax.set_xlim(-0.5, 1.75)
for t, col in zip(ax.get_xticklabels(), (EU, KC)): t.set_color(col)
fs.letter(ax, "e")
ax.set_ylabel("representatives meeting the rule\n(≥ 50% of an egg mass's frogs, ≥ 3)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=SH[k], edgecolor="black", linewidth=0.5, label="%s egg mass%s" % (k, "" if k == "1" else "es")) for k in ("3", "2", "1")],
          frameon=False, fontsize=5, loc="upper left", bbox_to_anchor=(1.0, 1.0), title="WF22 months 2-3", title_fontsize=5)
fs.save(fig, "fig3_wf22_eggmass_breadth")
