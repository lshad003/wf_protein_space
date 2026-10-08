# Step 25e: README.md, second n88 pass: Steps 3, 4, 5, 6, 13 (criteria), 16
# (mapping summaries), 17 (aggregate test) and 19 (complete-only lengths, which
# are unchanged). Guarded literal replacements: each old text must occur exactly
# once, checked before any write and again as each edit is applied.
# New numbers come from results/n88_number_table.tsv; 2,396 is from logs/step2d_n88.log.
# usage: step25e_readme_n88_part2.py [--apply]
import sys
P = "/bigdata/stajichlab/lshad003/wf_protein_space/README.md"
t = open(P).read()
E = []
def e(old, new): E.append((old, new))
def sup(path, purpose, tick, n88, n88purpose, extra=""):
    q = "`" if tick else ""
    n88p = path.replace(path.rsplit("/", 1)[1], n88)
    e("| %s%s%s | %s |" % (q, path, q, purpose),
      "| %s%s%s | %s; 89 metagenomes, superseded by %s%s%s |\n| %s%s%s | %s |" % (q, path, q, purpose, q, n88p, q, q, n88p, q, n88purpose) + extra)

# ---- Step 3
e("mapped proteins had median length 261 aa against 123 aa for unmapped", "mapped proteins had median length 260 aa against 123 aa for unmapped")
e("and on that setting 46.15% of new-cohort proteins match the earlier catalog at gene level and 66.66% at family level.",
  "and on that setting 45.89% of new-cohort proteins match the earlier catalog at gene level and 66.49% at family level. At 88 metagenomes the 2,396 UHM586.41010 proteins are dropped from the new-cohort subsample and the existing hits are restricted to the remaining 196,446; no search is repeated, and the WF22 control is unaffected.")
sup("scripts/step2d_control_and_tiers.sh", "Positive control and family-level tier", 1, "step2d_control_and_tiers_n88.sh", "Control and family-level tier, UHM586.41010 queries dropped from the existing hits")
sup("scripts/step2e_covmode2.sh", "Coverage-mode comparison with control", 1, "step2e_covmode2_n88.sh", "Coverage-mode comparison, UHM586.41010 queries dropped from the existing hits")

# ---- Step 4
e("**Main result.** 55.27% of gene-level clusters are singletons.", "**Main result.** 55.29% of gene-level clusters are singletons.")
e("clusters containing all three cohorts go from 6.71% to 20.06%", "clusters containing all three cohorts go from 6.72% to 20.06%")
e("top 1% of family-level clusters holding 45.85% of all proteins.", "top 1% of family-level clusters holding 45.77% of all proteins.")
sup("scripts/step3_cluster_composition.sh", "Cluster membership tables and cohort composition", 1, "step3_cluster_composition_n88.sh", "Cohort composition from the 88-metagenome cluster tables")
sup("scripts/step4_nonsingleton.sh", "Composition under four support filters, size distribution", 1, "step4_nonsingleton_n88.sh", "Composition under four support filters, 88 metagenomes; reproduces 2,067,011 supported clusters")

# ---- Step 5
e("WF22-exclusive clusters fall from 15.12% to 9.26% once the sample-count advantage is removed, while WF24-exclusive remains 26.13% and clusters shared by all three cohorts are 24.22%. Per-cohort richness at 9 metagenomes is 765,122 clusters for WF22, 1,166,662 for WF23 and 2,331,079 for WF24, tracking per-sample protein yield rather than cohort identity.",
  "WF22-exclusive clusters fall from 15.27% to 9.52% once the sample-count advantage is removed, while WF24-exclusive remains 25.28% and clusters shared by all three cohorts are 24.33%. Per-cohort richness at 9 metagenomes is 765,122 clusters for WF22, 1,166,662 for WF23 and 2,080,538 for WF24, tracking per-sample protein yield rather than cohort identity. Draws are seeded; with 35 WF24 metagenomes the same seeds select different WF24 animals, so the WF24 figures reflect a different draw as well as the removal of UHM586.41010.")
sup("scripts/step5_equaln_rarefaction.sh", "Equal-n cohort comparison and per-cohort accumulation", 1, "step5_equaln_rarefaction_n88.sh", "Equal-n comparison and per-cohort accumulation, 88 metagenomes, same seeds")

# ---- Step 6
e("Sequencing run, centre, platform and read depth are tabulated for all 89\nmetagenomes,", "Sequencing run, centre, platform and read depth are tabulated for all 88\nmetagenomes,")
e("Within the shared run WF23 and WF24 yield 422,857 and 415,208 proteins per metagenome", "Within the shared run WF23 and WF24 yield 422,857 and 431,504 proteins per metagenome")
e("Across its two runs WF24 yields 415,208 against 756,692 proteins per metagenome, a 1.8-fold difference on 1.24-fold more reads.",
  "Across its two runs WF24 yields 431,504 against 756,692 proteins per metagenome, a 1.75-fold difference on 1.24-fold more reads.")
sup("scripts/step6b_batch_table.sh", "Batch table and the two key contrasts", 1, "step6b_batch_table_n88.sh", "Batch table and the two key contrasts, 88 metagenomes")
e("Output: `results/batch_table_89.tsv`", "Output: `results/batch_table_89.tsv`; at 88: `results/batch_table_88.tsv`")

# ---- Step 13 criteria (the mobile-element sentence is left; see report)
e("runs 374,503 (63.88%) at E <= 1e-5 with no coverage requirement, 430,833\n(73.49%) at E <= 1e-10 with query and subject coverage at least 50%, and\n449,170 (76.62%) at E <= 1e-20",
  "runs 374,133 (63.92%) at E <= 1e-5 with no coverage requirement, 430,301\n(73.51%) at E <= 1e-10 with query and subject coverage at least 50%, and\n448,567 (76.63%) at E <= 1e-20")
sup("scripts/step13h_gu_eu_split.sh", "Split reported under several stated criteria", 0, "step13h_gu_eu_split_n88.sh", "Split under the same criteria, supported unknowns at 88 metagenomes",
    "\n| scripts/step13h_gu_eu_split_89log.sh | Unchanged step13h run once more to record its 89-sample output in a log |")

# ---- Step 16 (no remapping)
e("0.9634 in WF23 (n = 9) and 0.9671 in WF24 (n = 36), a spread under", "0.9634 in WF23 (n = 9) and 0.9667 in WF24 (n = 35), a spread under")
e("Genes detected at least once:\n2,069,187 under primary counting and 2,042,647 under the MAPQ >= 10 subset.",
  "Genes detected at least once:\n2,066,746 under primary counting and 2,040,274 under the MAPQ >= 10 subset.")
e("Reads were not remapped, so\nthe reference, mapping and verification figures above describe the 89-sample\nmapping.",
  "Reads were not remapped, so\nthe reference and the column-sum verification above describe the 89-sample\nmapping; the mapped fractions and detected-gene counts are given for the 88\nmetagenomes (results/mapping_qc_88.tsv, results/count_matrix_mapq10_n88.tsv.gz).")
sup("scripts/step16t_qc_table.sh", "Per-sample mapping quality table with cohort labels", 0, "step16t_qc_table_n88.sh", "Mapping quality table without UHM586.41010, by cohort")
sup("scripts/step16u_matrix.py", "Count matrices assembled and checked against the summaries", 0, "step16u_matrix_n88.py", "Both matrices at 88 metagenomes from the 89 matrices; column sums checked")

# ---- Step 17 aggregate test (WF24 at 34 animals)
e("differ by treatment in any cohort (Kruskal-Wallis p = 0.46 and 0.42 in WF22,\nMann-Whitney p = 0.11 and 0.19 in WF23, Kruskal-Wallis p = 0.81 and 0.85 in\nWF24).",
  "differ by treatment in any cohort (Kruskal-Wallis p = 0.4688 and 0.4222 in WF22,\nMann-Whitney p = 0.1111 and 0.1905 in WF23, Kruskal-Wallis p = 0.4481 and 0.5527\nin WF24 with 34 animals).")
sup("scripts/step17b_class_shift.py", "Unannotated abundance share tested by treatment within cohort", 0, "step17b_class_shift_n88.py", "Unannotated abundance share by treatment, 88 metagenomes, WF24 34 animals",
    "\n| scripts/step17b_class_shift_89log.sh | Unchanged step17b run once more to record its 89-sample output in a log |")

# ---- Step 19 (values unchanged; script table only)
sup("scripts/step16y_complete_by_class.py", "Reading frame completeness by class from the gene caller flags", 0, "step16y_complete_by_class_n88.py", "Completeness and mean length by class, 88 metagenomes")
e("results/core_unknown_by_years_n88.tsv.gz", "results/core_unknown_by_years_n88.tsv.gz, results/completeness_by_class_n88.tsv")

fail = [(t.count(o), o) for o, _ in E if t.count(o) != 1]
print("edits:", len(E), " failing guard:", len(fail))
for c, o in fail: print("  FOUND %d TIMES: %r" % (c, o[:110]))
if fail: sys.exit("REFUSED: README.md not written")
for o, n in E:
    if t.count(o) != 1: sys.exit("REFUSED: earlier edit changed the count of %r; README.md not written" % o[:110])
    t = t.replace(o, n)
if "--apply" in sys.argv:
    open(P, "w").write(t); print("wrote", P)
else:
    print("dry run, README.md not written")
