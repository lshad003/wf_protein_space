# Figure 3 additions (Steps 34t, 34t2): tank panel f, grouped bars of tank-restricted EU and K reps with the 95th percentile
# of each null, WF24 egg mass 4 and WF23 egg mass 3, assembly presence and >= 5 reads at c_id95 (plus a slide version of
# WF24 EM4 assembly); and a stacked bar of WF22 months 2-3 reps meeting the Step 34c3 rule in exactly 1, 2 or
# 3 egg masses (assembly; panel e).
# Inputs: results/step34t/tank_run_check.tsv, tank_run_crosstab.tsv, wf22_eggmass_breadth.tsv.
# Outputs: figures/fig3_tank_panel.*, figures/fig3_tank_panel_slide.*, figures/fig3_wf22_eggmass_breadth.*
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
fs.setup(); R = fs.W + "/results/step34t/"
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
EU, KC = "#AA3377", "#4477AA"
# Panel f (Step 34t2): grouped bars of observed tank-restricted reps; dark tick = 95th percentile of that class's null;
#   symbol = that class vs its own null; EU-K bracket only where the EU - K test has p <= 0.05. Final p-values from
#   results/step34t/tank_run_check.tsv: WF24 EM4 within-run permutation (1,000), WF23 EM3 exact (84 labelings). Refuses to
#   draw if results/step34t/tank_run_crosstab.tsv marks a contrast's tank and run as nested.
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
X = rows("tank_run_crosstab.tsv")
if any(r["tank"] == "nested" and r["n_froglets"] == "1" for r in X): raise SystemExit("REFUSING: tank and run nested")
T = rows("tank_run_check.tsv"); FINAL = {"WF24_EM4_tank": "within_run", "WF23_EM3_tank": "exact"}
def get(c, m, s): return next(r for r in T if r["contrast"] == c and r["presence"] == m and r["set"] == s and r["method"] == FINAL[c])
def star(p): return "***" if p <= 0.001 else "**" if p <= 0.01 else "*" if p <= 0.05 else "n.s."
NFR = {"WF24_EM4_tank": "13 vs 19", "WF23_EM3_tank": "3 vs 6"}
CON = [("WF24_EM4_tank", "WF24 EM4"), ("WF23_EM3_tank", "WF23 EM3")]
MODES = [("assembly", "assembly presence"), ("reads>=5_c_id95", "\u2265 5 reads (\u2265 95% identity)")]
def tank_panel(axes, cons, modes, fsz, bw=0.34):
    top = max(max(int(get(c, m, s)["observed_restricted"]), float(get(c, m, s)["null_p95"])) for c, _ in cons for m, _ in modes for s in ("EU", "K"))
    for ax, (mode, mlab) in zip(axes, modes):
        for gi, (c, lab) in enumerate(cons):
            hs = {}
            for k, (s, col) in enumerate((("EU", EU), ("K", KC))):
                r = get(c, mode, s); x = gi + (-0.5 + k) * bw; o = int(r["observed_restricted"]); p95 = float(r["null_p95"])
                ax.bar(x, o, bw * 0.92, color=col, linewidth=0)
                ax.hlines(p95, x - bw * 0.46, x + bw * 0.46, color="#222222", linewidth=1.4, zorder=3)
                h = max(o, p95) + 0.02 * top; hs[s] = h
                ax.text(x, h, "%s\n%s" % (format(o, ","), star(float(r["p_one_sided"]))), ha="center", va="bottom", fontsize=fsz, linespacing=1.1)
            d = get(c, mode, "D_EU_minus_K")
            if float(d["p_one_sided"]) <= 0.05:
                y = max(hs.values()) + 0.17 * top; x0, x1 = gi - 0.5 * bw, gi + 0.5 * bw
                ax.plot([x0, x0, x1, x1], [y - 0.025 * top, y, y, y - 0.025 * top], color="black", linewidth=0.6)
                ax.text(gi, y, star(float(d["p_one_sided"])), ha="center", va="bottom", fontsize=fsz)
        ax.set_xticks(range(len(cons))); ax.set_xticklabels(["%s\n(%s froglets)" % (lab, NFR[c]) for c, lab in cons], fontsize=fsz)
        ax.set_xlim(-0.6, len(cons) - 0.4); ax.set_ylim(0, top * 1.38); ax.set_title(mlab, fontsize=fsz + 1)
    axes[0].set_ylabel("tank-restricted representatives", fontsize=fsz + 0.5)
    return [Patch(color=EU, label="EU"), Patch(color=KC, label="K"), Line2D([], [], color="#222222", linewidth=1.4, label="chance, 95th percentile")]
fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.2), sharey=True, gridspec_kw=dict(wspace=0.12))
hl = tank_panel(axes, CON, MODES, 5)
axes[1].tick_params(axis="y", labelleft=False)
fig.legend(handles=hl, frameon=False, fontsize=5, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.08))
fs.letter(axes[0], "f")
fs.save(fig, "fig3_tank_panel")
plt.rcParams.update({"font.size": 11, "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.linewidth": 0.8})
fig, ax = plt.subplots(figsize=(3.6, 3.8))
hl = tank_panel([ax], CON[:1], MODES[:1], 11, bw=0.36)
ax.set_xlim(-0.5, 0.5); ax.legend(handles=hl, frameon=False, fontsize=9, loc="upper center", ncol=1, bbox_to_anchor=(0.5, -0.18))
fs.save(fig, "fig3_tank_panel_slide")
fs.setup()
print("\nLEGEND DRAFT, panel f")
print("Tank restriction: a representative is tank-restricted when present in >= 50% of one tank's froglets and in >= 3 of them, and")
print("absent from every froglet of the other tank of the same egg mass; counts are summed over both tanks. Presence = any assembled")
print("member (assembly) or >= 5 reads at >= 95% identity (c_id95). 41,114 EU and 41,114 randomly drawn K representatives (seed 343).")
print("Chance = tank labels reassigned among the egg mass's froglets with tank sizes kept: WF24 EM4, 1,000 permutations within")
print("sequencing run; WF23 EM3, all 84 possible labelings (exact). p one-sided; *** <= 0.001, ** <= 0.01, * <= 0.05, n.s. > 0.05.")
print("Bracket: EU - K difference vs the same permutations, shown only where p <= 0.05. Dark tick: 95th percentile of the class's null.")
for c, lab in CON:
    for m, ml in MODES:
        for s in ("EU", "K", "D_EU_minus_K"):
            r = get(c, m, s)
            print("  %-9s %-28s %-13s obs %5s  p %.4f (%s)  null mean %s  sd %s  p95 %s  [%s, n=%s]" % (lab, ml, s, r["observed_restricted"], float(r["p_one_sided"]),
                  star(float(r["p_one_sided"])), r["null_mean"], r["null_sd"], r["null_p95"], r["method"], r["n_null"]))
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
ax.set_ylabel("representatives meeting the rule\n(≥ 50% of an egg mass's froglets, ≥ 3)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=SH[k], edgecolor="black", linewidth=0.5, label="%s egg mass%s" % (k, "" if k == "1" else "es")) for k in ("3", "2", "1")],
          frameon=False, fontsize=5, loc="upper left", bbox_to_anchor=(1.0, 1.0), title="WF22 months 2-3", title_fontsize=5)
fs.save(fig, "fig3_wf22_eggmass_breadth")
