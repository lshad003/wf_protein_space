# Step 25d: README.md from 89 to 88 metagenomes (prompts/n88_readme.md).
# Every edit is a guarded literal replacement: the old text must occur exactly
# once, otherwise nothing is written and the failures are reported.
# New numbers come from results/n88_number_table.tsv; the few that do not are
# marked SRC in the comment beside the edit.
# usage: step25d_readme_n88.py [--apply]
import sys
W = "/bigdata/stajichlab/lshad003/wf_protein_space/"
P = W + "README.md"
t = open(P).read()
E = []
def e(old, new): E.append((old, new))
def sup(path, purpose, tick, n88, n88purpose):
    q = "`" if tick else ""
    n88p = path.replace(path.rsplit("/", 1)[1], n88)
    e("| %s%s%s | %s |" % (q, path, q, purpose),
      "| %s%s%s | %s; 89 metagenomes, superseded by %s%s%s |\n| %s%s%s | %s |" % (q, path, q, purpose, q, n88p, q, q, n88p, q, n88purpose))

# ---- Sample set (rule 2). SRC: metagenome and animal counts from metadata/samples_88_treatment.tsv
# and metadata/wf22_design.tsv; 9,637 from logs/step25b.log.
e("89 faecal metagenomes from 64 animals. Faecal", "88 faecal metagenomes from 59 animals. Faecal")
e("| WF24 | 36 | 36 | one monthly pool per animal |", "| WF24 | 35 | 35 | one monthly pool per animal |")
e("""month. Tank is confounded
with egg mass in WF24.
""", """month. Tank is confounded
with egg mass in WF24.

**Catalog at 88 metagenomes.** UHM586.41010 (WF24 code 4 control, elevated
*Basidiobolus* reads per J. Dallas) was removed after clustering. The MMseqs2
clustering was run on all 89 metagenomes and is not repeated: every member
protein from UHM586.41010 is removed from the cluster tables, the support rule
(at least three members from at least two metagenomes) is reapplied, and the
original representatives are kept, including 9,637 clusters at 95% whose
representative sequence comes from UHM586.41010. Reads are not remapped; the
UHM586.41010 column is dropped from the count matrix. The exclusion criterion
is the one applied to WF22 control UHM56.10839, which was excluded before the
catalog was built. UHM590.41012 stays in the catalog and is excluded only from
the treatment tests, so WF24 treatment tests use 34 animals. Every number
changed by the removal is listed with its old value and output file in
results/n88_number_table.tsv. The sequenced WF24 description above (36
animals) is unchanged.

| File | Purpose |
|---|---|
| scripts/step25a_filter_n88.py | Supported set, count matrix and sample table at 88 metagenomes |
| scripts/step25b_membership_n88.py | Cluster tables at 88 metagenomes, all three tiers; unknown families restricted to supported members |
| scripts/step25c_number_table.py | Old and new value of every README number, parsed from logs and outputs |
| scripts/step24b_accum_n88.py | Gene accumulation against samples at 88 metagenomes (89 version step24b_accum.py, no README step) |

Output: results/supported_reps_95_n88.txt, results/count_matrix_primary_n88.tsv.gz,
metadata/samples_88_treatment.tsv, results/clusters_{95,90,50}_n88.tsv,
results/reps_95_n88.txt, results/unk_clusters_50_n88.tsv,
results/accumulation_by_sample_n88.tsv, results/n88_number_table.tsv
""")

# ---- Step 1 (rule 3). SRC: 292,633 from results/step1_protein_counts_45.tsv
e("UHM585.41009 is a\nlow outlier at 313,220 proteins from a 64 Mb assembly.",
  "UHM585.41009 is a\nlow outlier at 313,220 proteins from a 64 Mb assembly. UHM586.41010 (292,633\nproteins, the lowest WF24 yield) was removed from the catalog afterwards; see\nSample set.")

# ---- Step 2: clustering stays at 89; state what remains (table rows: pooled clusters, 50% families all)
e("2,922,537 at 50%. Catalog is `catalog/db/LsPS_AA`.",
  "2,922,537 at 50%. Catalog is `catalog/db/LsPS_AA`. After removing the\nUHM586.41010 members, 6,171,602 clusters remain at 95% and 2,918,630 at 50%.")

# ---- Step 4 (table: supported total)
e("metagenomes leaves 2,069,453 clusters", "metagenomes leaves 2,067,011 clusters")

# ---- Step 7
e("WF24 gives 1,285,504 in one run", "WF24 gives 1,295,364 in one run")
e("differs by 31% across runs while different cohorts within one run differ by 11%.",
  "differs by 30.3% across runs while different cohorts within one run differ by 11.7%.")
e("reaches 6,182,117 clusters at 40,550,595 proteins and is still adding roughly 435,000 clusters per 5 million proteins",
  "reaches 6,171,602 clusters at 40,257,962 proteins and is still adding 436,425 clusters between 35 and 40 million proteins")
sup("scripts/step7_rarefaction.sh", "Accumulation against proteins sampled, per cohort and pooled", 1, "step7_rarefaction_n88.sh", "Accumulation against proteins sampled, 88 metagenomes")
sup("scripts/step7b_wf24_byrun.sh", "Richness at matched depth, WF24 split by run", 1, "step7b_wf24_byrun_n88.sh", "Richness at matched depth, WF24 split by run, 88 metagenomes")
e("Output: `results/rarefaction_95.tsv`", "Output: `results/rarefaction_95.tsv`; at 88: `results/rarefaction_95_n88.tsv`")

# ---- Step 9
e("46.91% of the 40,550,595 predicted proteins", "46.89% of the 40,257,962 predicted proteins")
e("is lower at 35.28%", "is lower at 35.33%")
e("11.48% for singletons, 33.13% for clusters of two, 61.07% for three to five, 81.94% for six to twenty, and 86.14% for twenty-one",
  "11.58% for singletons, 32.98% for clusters of two, 61.16% for three to five, 81.99% for six to twenty, and 86.11% for twenty-one")
e("support-filtered set of 2,069,453 clusters as the analytical catalog", "support-filtered set of 2,067,011 clusters as the analytical catalog")
sup("scripts/step9_orf_completeness.sh", "Partial flags extracted, completeness against cluster size", 1, "step9_orf_completeness_n88.sh", "Completeness against cluster size, 88 metagenomes")
e("Output: `results/orf_completeness.tsv`", "Output: `results/orf_completeness.tsv`; at 88: `results/orf_completeness_n88.tsv`")

# ---- Step 10
e("The 6,182,117 gene-level representatives are annotated", "The 6,171,602 gene-level representatives are annotated")
e("**Main result.** 48.27% of gene-level", "**Main result.** 48.29% of gene-level")
e("from 36.79% for singletons to 90.51% for", "from 36.83% for singletons to 90.53% for")
e("(65.1% against 40.1% among singletons)", "(65.1% against 40.4% among singletons)")
e("76.01% of all families", "76.03% of all families")
e("falls to 55.35% of 1,022,162 families", "falls to 55.40% of 1,020,958 families")
e("leaving 565,769 supported", "leaving 565,611 supported")
sup("scripts/step10_annotation_summary.sh", "Annotated fraction by cluster size and completeness", 1, "step10_annotation_summary_n88.sh", "Annotated fraction by cluster size and completeness, 88 metagenomes")
sup("scripts/step11_interaction_and_supported.sh", "Completeness within cluster size, family-tier dark fraction under support filters", 1, "step11_interaction_and_supported_n88.sh", "Family-tier dark fraction under support filters, 88 metagenomes")
e("Output: `results/annotation_summary.tsv`, `results/family_dark_fraction.tsv`",
  "Output: `results/annotation_summary.tsv`, `results/family_dark_fraction.tsv`; at 88:\n`results/annotation_summary_n88.tsv`, `results/family_dark_fraction_n88.tsv`")

# ---- Step 11
e("0.16% Viruses of 2,984,013 annotated", "0.16% Viruses of 2,980,052 annotated")
e("Bacteroidetes 13.94%", "Bacteroidetes 13.95%")
e("Metazoa at 3.04%", "Metazoa at 3.05%")
e("from 7.64% among singletons to 3.59% among", "from 7.63% among singletons to 3.61% among")
sup("scripts/step13_taxonomy.sh", "Broad clade and assignment level, against cluster size", 1, "step13_taxonomy_n88.sh", "Broad clade and assignment level, 88 metagenomes")
e("Output: `results/taxonomy_summary.tsv`", "Output: `results/taxonomy_summary.tsv`; at 88: `results/taxonomy_summary_n88.tsv`")

# ---- Step 12
e("""**Coverage.** Of 6,182,117 representatives, eggNOG annotates 2,984,013
(48.27%), Pfam 2,713,156 (43.89%), and their union 3,227,284 (52.20%). Within
the support-filtered set of 2,069,453, the union leaves 586,215 unannotated
(28.33%). At the 50% family tier the supported dark fraction falls from 55.35%
under eggNOG alone to 51.42% under the union""",
"""**Coverage.** Of 6,171,602 representatives, eggNOG annotates 2,980,052
(48.29%), Pfam 2,709,544 (43.90%), and their union 3,222,677 (52.22%). Within
the support-filtered set of 2,067,011, the union leaves 585,347 unannotated
(28.32%). At the 50% family tier the supported dark fraction falls from 55.40%
under eggNOG alone to 51.47% under the union""")
e("""counts are K 1,344,403, KWP 138,835, U 586,215. Of the K representatives,
76,589 carry only DUF or UPF domains, which a stricter definition would move
out of K, leaving 1,267,814.""",
"""counts are K 1,343,020, KWP 138,644, U 585,347. Of the K representatives,
76,515 carry only DUF or UPF domains, which a stricter definition would move
out of K, leaving 1,266,505.""")
sup("scripts/step12d_union_replevel.sh", "Representative-level union of eggNOG and Pfam", 0, "step12d_union_replevel_n88.sh", "Representative-level union, 88 metagenomes")
sup("scripts/step12e_family_union.py", "Union carried to the 50% family tier", 0, "step12e_family_union_n88.py", "Union carried to the 50% family tier, 88 metagenomes")
sup("scripts/step12f_rep_classes.sh", "Per-representative K, KWP and U labels with the DUF split", 0, "step12f_rep_classes_n88.sh", "K, KWP and U labels, 88 metagenomes; classes checked equal to the 89 labels")
e("""Output: results/union_rep_level.tsv, results/family_union_dark.tsv,
results/rep_classes.tsv, results/known_union_ids.txt""",
"""Output: results/union_rep_level.tsv, results/family_union_dark.tsv,
results/rep_classes.tsv, results/known_union_ids.txt; at 88:
results/union_rep_level_n88.tsv, results/family_union_dark_n88.tsv,
results/rep_classes_n88.tsv, results/known_union_ids_n88.txt""")

# ---- Step 13 (search itself ran on the 89 set and is kept as described)
e("Thirteen hours on 32 cores with the database staged to node-local scratch.",
  "Thirteen hours on 32 cores with the database staged to node-local scratch.\nAt 88 metagenomes 585,347 of these remain supported and are classified from\nthe same search.")
e("Of the supported catalogue, 430,831 representatives (20.82%)", "Of the supported catalogue, 430,299 representatives (20.82%)")
e("155,384 representatives (7.51%)", "155,048 representatives (7.50%)")
sup("scripts/step13i_fourway.py", "Four-way classification assembled with per-class properties", 0, "step13i_fourway_n88.py", "Four-way classification, 88 metagenomes")
e("Output: results/step13_nrclust20260128_hits.tsv, results/fourway_classes.tsv",
  "Output: results/step13_nrclust20260128_hits.tsv, results/fourway_classes.tsv; at 88:\nresults/fourway_classes_n88.tsv")

# ---- Step 14 (known-gene numbers and the 40% control left as they are, rule 7)
e("AntiFam v6.0 flags 1,477 of 586,215\nrepresentatives, 0.2520%.", "AntiFam v6.0 flags 1,464 of 585,347\nrepresentatives, 0.2501%.")
e("coverage, the 586,215 unknown representatives give 536,813 families. Only 8,482\nhave three or more members and 10 have 100 or more.",
  "coverage, the 585,347 unknown representatives give 536,030 families. Only 8,464\nhave three or more members and 10 have 100 or more. At 88 metagenomes these are\nthe families of the 89-sample clustering restricted to supported members, not a\nnew clustering.")
e("cannot be\nbuilt from 89.", "cannot be\nbuilt from 88.")
sup("scripts/step14b_antifam.sh", "AntiFam screen of the unannotated representatives", 0, "step14b_antifam_n88.sh", "AntiFam screen rerun on the 585,347 supported unknowns at 88 metagenomes")
sup("scripts/step14c_novel_families.sh", "Unknown representatives clustered at 50% identity", 0, "step14c_novel_families_n88.sh", "Family sizes at 88 metagenomes from the 50% clustering, no reclustering")
e("""Output: results/antifam_hits.tblout, results/unk_clusters_50.tsv,
results/unk_family_sizes.tsv""",
"""Output: results/antifam_hits.tblout, results/unk_clusters_50.tsv,
results/unk_family_sizes.tsv; at 88: results/antifam_hits_n88.tblout,
results/unk_clusters_50_n88.tsv, results/unk_family_sizes_n88.tsv""")

# ---- Step 16 (mapping not repeated)
e("""Output: results/count_matrix_primary.tsv.gz, results/count_matrix_mapq10.tsv.gz,
results/mapping_qc_89.tsv, results/supported_reps_95.txt""",
"""Output: results/count_matrix_primary.tsv.gz, results/count_matrix_mapq10.tsv.gz,
results/mapping_qc_89.tsv, results/supported_reps_95.txt

At 88 metagenomes the analyses use results/count_matrix_primary_n88.tsv.gz:
the UHM586.41010 column removed and rows restricted to the 2,067,011 clusters
still supported (scripts/step25a_filter_n88.py). Reads were not remapped, so
the reference, mapping and verification figures above describe the 89-sample
mapping.""")

# ---- Step 17 (LinDA sentence keeps 368,236: that is the size LinDA was run on)
e("contrasts return 4,993 and 7,275 genes.", "contrasts return 5,005 and 7,263 genes.")
e("returns 57,183 and 51,864 genes and month returns up to 141,713.", "returns 57,120 and 51,808 genes and month returns up to 141,580.")
sup("scripts/step17g_limma_wf22.R", "WF22 tested blocked and unblocked, and within egg mass 3", 0, "step17g_limma_wf22_n88.R", "WF22 tested blocked and unblocked, and within egg mass 3, 88-metagenome matrix")
sup("scripts/step17h_limma_wf2324.R", "WF23 treatment contrasts; its WF24 run (36 animals) is superseded by step17i", 0, "step17h_limma_wf2324_n88.R", "WF23 treatment contrasts, 88-metagenome matrix")
sup("scripts/step17i_limma_wf24_n34.R", "WF24 treatment contrasts, 34 animals, exclusions in metadata/wf24_excluded.tsv", 0, "step17i_limma_wf24_n34_n88.R", "WF24 treatment contrasts, 34 animals, 88-metagenome matrix")
e("metadata/wf22_design.tsv, metadata/wf24_treatment_key.tsv, metadata/wf24_excluded.tsv",
  "metadata/wf22_design.tsv, metadata/wf24_treatment_key.tsv, metadata/wf24_excluded.tsv;\nat 88: results/limma_*_n88.csv, metadata/samples_88_treatment.tsv")

# ---- Step 19
e("795 for KWP, 497 for GU", "795 for KWP, 498 for GU")
e("up 79.5% of K, 69.7% of KWP, 75.9% of GU", "up 79.5% of K, 69.8% of KWP, 76.0% of GU")
e("**They are widespread.** 82.9% of unannotated genes occur in more than 20 of\n89 samples and 3,471 occur in all 89. Only 0.7% occur in five or fewer.",
  "**They are widespread.** 82.77% of unannotated genes occur in more than 20 of\n88 samples and 3,496 occur in all 88. Only 0.65% occur in five or fewer.")
e("""497,044 unannotated genes occur in all three
collection years and carry 91.18% of unannotated gene reads. The 18,914
single-year genes are shorter, mean 246 bp, and carry 1.22%. A further 261
genes""", """496,266 unannotated genes occur in all three
collection years and carry 91.15% of unannotated gene reads. The 18,896
single-year genes are shorter, mean 246 bp, and carry 1.23%. A further 260
genes""")
e("Unannotated genes are 28.33% of\nsupported genes but 7.03% of raw mapped reads. After length normalization they\nare 19.28% of gene copies.",
  "Unannotated genes are 28.32% of\nsupported genes but 7.04% of raw mapped reads. After length normalization they\nare 19.32% of gene copies.")
e("7.11%\nin WF22, 6.29% in WF23 and 7.11% in WF24, with a per-sample range of 4.39 to\n16.13%.",
  "7.11%\nin WF22, 6.28% in WF23 and 7.16% in WF24, with a per-sample range of 4.38 to\n16.12%.")
sup("scripts/step16v_prevalence.py", "Prevalence per gene by cohort and class", 0, "step16v_prevalence_n88.py", "Prevalence per gene by cohort and class, 88 metagenomes")
sup("scripts/step16w_dark_abundance.py", "Abundance share of each class, raw and length-normalized", 0, "step16w_dark_abundance_n88.py", "Abundance share of each class, 88 metagenomes")
sup("scripts/step19a_core_unknown.py", "Unannotated genes characterized by number of years detected", 0, "step19a_core_unknown_n88.py", "Unannotated genes by number of years detected, 88 metagenomes")
e("figures fig_prevalence_by_class.png and fig_dark_abundance.png",
  "figures fig_prevalence_by_class.png and fig_dark_abundance.png; at 88:\nresults/prevalence_primary_n88.tsv.gz, results/dark_abundance_by_sample_n88.tsv,\nresults/core_unknown_by_years_n88.tsv.gz")

# ---- Step 20
e("Contig assignments exist for 48 of 89 samples, giving taxonomy\nfor 888,392 of 2,069,453 genes, about 43%.",
  "Contig assignments exist for 48 of 88 samples, giving taxonomy\nfor 888,183 of 2,067,011 genes, 42.97%.")
e("74.4% of GU and 91.5% of EU", "74.4% of GU and 91.6% of EU")
sup("scripts/step20c_gene_tax.py", "Contig taxonomy inherited by gene and summarized by class", 0, "step20c_gene_tax_n88.py", "Contig taxonomy by gene, 88 metagenomes; reads .tsv or .tsv.gz contig files")
e("Output: results/gene_taxonomy.tsv", "Output: results/gene_taxonomy.tsv; at 88: results/gene_taxonomy_n88.tsv")

# ---- Step 21 (rule 4)
e("applied only to the 8,482 unknown families", "applied only to the 8,464 unknown families")
e("Of the 8,482 testable families,\n3,628 (42.8%)", "Of the 8,464 testable families,\n3,614 (42.7%)")
e("which\nis 1.4% of unannotated genes,", "whose\nmembers are 6.29% of supported unannotated genes,")
sup("scripts/step21a_neighbour_check.sh", "Gene identifiers confirmed to encode contig and position", 0, "step21a_neighbour_check_n88.sh", "Identifier check on the 88-metagenome tables")
sup("scripts/step21b_neighbours.py", "Annotated domains counted beside members of each unknown family", 0, "step21b_neighbours_n88.py", "Neighbour domains of unknown families, 88 metagenomes")
e("Output: results/unk_family_neighbours.tsv", "Output: results/unk_family_neighbours.tsv; at 88: results/unk_family_neighbours_n88.tsv")

# ---- Step 22 (R2 percentages and unchanged p values left as written)
e("treatment R2 is 0.173 (p = 0.678)", "treatment R2 is 0.173 (p = 0.679)")
e("gives 1,639,294 genes that vary across samples, a mean of 523,670 per sample.",
  "gives 1,637,544 genes that vary across samples, a mean of 523,226 per sample.")
e("does not differ by treatment (p = 0.356).", "does not differ by treatment (p = 0.3676).")
e("59,297 representatives, of which 7,902 are present", "59,297 representatives, of which 7,891 are present")
e("1,763 and 2,587 for month", "1,761 and 2,586 for month")
e("abundance does not differ (p = 0.305).", "abundance does not differ (p = 0.3181).")
sup("scripts/step22b_permanova_fix.R", "Composition tested under three permutation designs", 0, "step22b_permanova_fix_n88.R", "Composition under three permutation designs, 88-metagenome matrix")
sup("scripts/step22c_jaccard.R", "Presence/absence composition and gene richness", 0, "step22c_jaccard_n88.R", "Presence/absence composition and gene richness, 88-metagenome matrix")
sup("scripts/step22e_permanova_wf24_n34.R", "WF24 abundance and presence/absence composition, 34 animals", 0, "step22e_permanova_wf24_n34_n88.R", "WF24 composition, 34 animals, 88-metagenome matrix")
sup("scripts/step22d_cazy_subset.R", "Carbohydrate subset tested per gene and by composition", 0, "step22d_cazy_subset_n88.R", "Carbohydrate subset per gene and by composition, 88-metagenome matrix")
e("Output: logs step22b, step22c, step22d and step22e; results/limma_cazy_*.csv",
  "Output: logs step22b, step22c, step22d and step22e; results/limma_cazy_*.csv;\nat 88: logs step22b_n88 to step22e_n88, results/limma_cazy_*_n88.csv")

fail = []
for old, new in E:
    c = t.count(old)
    if c != 1: fail.append((c, old))
print("edits:", len(E), " failing guard:", len(fail))
for c, old in fail: print("  FOUND %d TIMES: %r" % (c, old[:110]))
if fail: sys.exit("REFUSED: README.md not written")
for old, new in E:
    if t.count(old) != 1: sys.exit("REFUSED: earlier edit changed the count of %r; README.md not written" % old[:110])
    t = t.replace(old, new)
if "--apply" in sys.argv:
    open(P, "w").write(t); print("wrote", P)
else:
    print("dry run, README.md not written")
