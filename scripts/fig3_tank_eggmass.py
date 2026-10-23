# Figure 3 panels d, e, f (Steps 34t, 34t2, 34t5), as functions imported by scripts/fig3_recurrence.py (six-panel Fig. 3);
# run directly it writes the standalone panels.
# d: WF22 egg-mass-restricted EU and K reps per egg mass (assembly, months 2-3), bars with a dark tick at the 95th
#    percentile of that class's Step 34c3 null and one symbol for that class vs its own null; overall EU - K (Step 34c3b).
# e: stacked bar of WF22 months 2-3 reps meeting the Step 34c3 rule in exactly 1, 2 or 3 egg masses (assembly).
# f: tank panel, grouped bars of tank-restricted EU and K reps with the 95th percentile of each null, WF24 egg mass 4 and
#    WF23 egg mass 3, assembly presence and >= 5 reads at c_id95 (plus a slide version of WF24 EM4 assembly).
# Inputs: results/step34t/wf22_eggmass_nulls.tsv, tank_run_check.tsv, tank_run_crosstab.tsv, wf22_eggmass_breadth.tsv.
# Outputs (run directly): figures/fig3_tank_panel.*, figures/fig3_tank_panel_slide.*, figures/fig3_wf22_eggmass_breadth.*
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
R = fs.W + "/results/step34t/"
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
EU, KC = "#AA3377", "#4477AA"
def star(p): return "***" if p <= 0.001 else "**" if p <= 0.01 else "*" if p <= 0.05 else "n.s."
TICK = dict(color="#222222", linewidth=1.4, zorder=3)
def handles(): return [Patch(color=EU, label="EU"), Patch(color=KC, label="K"), Line2D([], [], label="chance, 95th percentile", **{k: v for k, v in TICK.items() if k != "zorder"})]
def bars(ax, x, o, p95, p, col, bw, top, fsz):
    ax.bar(x, o, bw * 0.92, color=col, linewidth=0)
    ax.hlines(p95, x - bw * 0.46, x + bw * 0.46, **TICK)
    h = max(o, p95) + 0.02 * top
    ax.text(x, h, "%s\n%s" % (format(o, ","), star(p)), ha="center", va="bottom", fontsize=fsz, linespacing=1.1)
    return h
# Panel f (Step 34t2): final p-values from results/step34t/tank_run_check.tsv: WF24 EM4 within-run permutation (1,000),
#   WF23 EM3 exact (84 labelings). Refuses to draw if results/step34t/tank_run_crosstab.tsv marks tank and run as nested.
X = rows("tank_run_crosstab.tsv")
if any(r["tank"] == "nested" and r["n_froglets"] == "1" for r in X): raise SystemExit("REFUSING: tank and run nested")
T = rows("tank_run_check.tsv"); FINAL = {"WF24_EM4_tank": "within_run", "WF23_EM3_tank": "exact"}
def get(c, m, s): return next(r for r in T if r["contrast"] == c and r["presence"] == m and r["set"] == s and r["method"] == FINAL[c])
NFR = {"WF24_EM4_tank": "13 vs 19", "WF23_EM3_tank": "3 vs 6"}
CON = [("WF24_EM4_tank", "WF24 EM4"), ("WF23_EM3_tank", "WF23 EM3")]
MODES = [("assembly", "assembly presence"), ("reads>=5_c_id95", "≥ 5 reads (≥ 95% identity)")]
def tank_panel(axes, cons, modes, fsz, bw=0.34):
    top = max(max(int(get(c, m, s)["observed_restricted"]), float(get(c, m, s)["null_p95"])) for c, _ in cons for m, _ in modes for s in ("EU", "K"))
    for ax, (mode, mlab) in zip(axes, modes):
        for gi, (c, lab) in enumerate(cons):
            hs = {}
            for k, (s, col) in enumerate((("EU", EU), ("K", KC))):
                r = get(c, mode, s)
                hs[s] = bars(ax, gi + (-0.5 + k) * bw, int(r["observed_restricted"]), float(r["null_p95"]), float(r["p_one_sided"]), col, bw, top, fsz)
            d = get(c, mode, "D_EU_minus_K")
            if float(d["p_one_sided"]) <= 0.05:
                y = max(hs.values()) + 0.26 * top; x0, x1 = gi - 0.5 * bw, gi + 0.5 * bw
                ax.plot([x0, x0, x1, x1], [y - 0.025 * top, y, y, y - 0.025 * top], color="black", linewidth=0.6)
                ax.text(gi, y, star(float(d["p_one_sided"])), ha="center", va="bottom", fontsize=fsz)
        ax.set_xticks(range(len(cons))); ax.set_xticklabels(["%s\n(%s froglets)" % (lab, NFR[c]) for c, lab in cons], fontsize=fsz)
        ax.tick_params(axis="y", labelsize=fsz)
        ax.set_xlim(-0.6, len(cons) - 0.4); ax.set_ylim(0, top * 1.5); ax.set_title(mlab, fontsize=fsz)
    axes[0].set_ylabel("tank-restricted representatives", fontsize=fsz)
    return handles()
# Panel d (Step 34t5): per-class nulls "step34c3" (as Step 34c3), overall EU - K from the shared permutations (34c3b).
N = rows("wf22_eggmass_nulls.tsv")
def getd(n, s, e): return next(r for r in N if r["null"] == n and r["set"] == s and r["egg_mass"] == e)
EMS = [("WF22_EM1", "EM1"), ("WF22_EM2", "EM2"), ("WF22_EM3", "EM3")]
def eggmass_panel(ax, fsz, bw=0.34):
    top = max(max(int(getd("step34c3", s, e)["observed_restricted"]), float(getd("step34c3", s, e)["null_p95"])) for e, _ in EMS for s in ("EU", "K"))
    for gi, (e, lab) in enumerate(EMS):
        for k, (s, col) in enumerate((("EU", EU), ("K", KC))):
            r = getd("step34c3", s, e)
            bars(ax, gi + (-0.5 + k) * bw, int(r["observed_restricted"]), float(r["null_p95"]), float(r["p_one_sided"]), col, bw, top, fsz)
    d = getd("shared", "D_EU_minus_K", "all")
    ax.text(0.02, 0.98, "EU − K, all egg masses: %s, p = %.3f" % (format(int(d["observed_restricted"]), ","), float(d["p_one_sided"])),
            transform=ax.transAxes, ha="left", va="top", fontsize=fsz)
    ax.set_xticks(range(len(EMS))); ax.set_xticklabels(["%s\n(%s froglets)" % (lab, getd("step34c3", "EU", e)["n_froglets"]) for e, lab in EMS], fontsize=fsz)
    ax.tick_params(axis="y", labelsize=fsz)
    ax.set_xlim(-0.6, len(EMS) - 0.4); ax.set_ylim(0, top * 1.3); ax.set_ylabel("egg-mass-restricted representatives", fontsize=fsz)
    return handles()
# Panel e
B = rows("wf22_eggmass_breadth.tsv")
SH = {"1": "#DDDDDD", "2": "#999999", "3": "#333333"}
def breadth_panel(ax, fsz):
    for i, (s, col) in enumerate((("EU", EU), ("K", KC))):
        bot = 0
        for k in ("1", "2", "3"):
            r = next(r for r in B if r["set"] == s and r["n_egg_masses_meeting_rule"] == k); n = int(r["n_reps"])
            ax.bar(i, n, 0.6, bottom=bot, color=SH[k], edgecolor=col, linewidth=1.0)
            if n: ax.text(i + 0.33, bot + n / 2, format(n, ","), va="center", fontsize=fsz)
            bot += n
    ax.set_xticks([0, 1]); ax.set_xticklabels(["EU", "K"], fontsize=fsz); ax.set_xlim(-0.5, 1.75); ax.tick_params(axis="y", labelsize=fsz)
    for t, col in zip(ax.get_xticklabels(), (EU, KC)): t.set_color(col)
    ax.set_ylabel("representatives meeting the rule\n(≥ 50% of an egg mass's froglets, ≥ 3)", fontsize=fsz)
    ax.legend(handles=[Patch(facecolor=SH[k], edgecolor="black", linewidth=0.5, label=k) for k in ("3", "2", "1")],
              frameon=False, fontsize=fsz, loc="upper left", bbox_to_anchor=(1.0, 1.0), title="egg masses\nmeeting the rule", title_fontsize=fsz)
def legend_d():
    print("\nLEGEND DRAFT, panel d")
    print("Egg-mass restriction, WF22 (15 froglets; months 2 and 3 pooled, presence = any assembled member in either sample): a")
    print("representative is egg-mass-restricted when present in >= 50% of one egg mass's froglets and in >= 3 of them, and absent from")
    print("every froglet of the other two WF22 egg masses. 41,114 EU and 41,114 randomly drawn K representatives (seed 343).")
    print("Chance = egg-mass labels reassigned among the 15 froglets with group sizes kept (3, 3, 9), 1,000 permutations per class.")
    print("Dark tick: 95th percentile of that class's null. Symbol: that class vs its own null, one-sided; *** <= 0.001, ** <= 0.01,")
    print("* <= 0.05, n.s. > 0.05. EU - K, all egg masses: the difference summed over egg masses, tested on 1,000 permutations shared by")
    print("EU and K, one-sided.")
    for e, lab in EMS + [("all", "all")]:
        for s in ("EU", "K"):
            r = getd("step34c3", s, e)
            print("  %-4s %-3s obs %4s  p %.4f (%s)  null mean %s  sd %s  p95 %s" % (lab, s, r["observed_restricted"], float(r["p_one_sided"]),
                  star(float(r["p_one_sided"])), r["null_mean"], r["null_sd"], r["null_p95"]))
    for e, lab in EMS + [("all", "all")]:
        r = getd("shared", "D_EU_minus_K", e)
        print("  %-4s EU - K %4s  p %.4f  null mean %s  sd %s  p95 %s (shared permutations)" % (lab, r["observed_restricted"], float(r["p_one_sided"]), r["null_mean"], r["null_sd"], r["null_p95"]))
def legend_f():
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
if __name__ == "__main__":
    fs.setup()
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
    fig, ax = plt.subplots(figsize=(1.9, 2.0))
    breadth_panel(ax, 5); fs.letter(ax, "e")
    fs.save(fig, "fig3_wf22_eggmass_breadth")
    legend_d(); legend_f()
