# Wood frog gut protein space: analysis record

This repository documents the construction of a non-redundant protein catalog
from wood frog (*Lithobates sylvaticus*) faecal metagenomes collected over three
years, and the analyses that characterize it: size and saturation, cohort
contributions, technical drivers of richness, and functional annotation.

## Sample set

88 faecal metagenomes from 59 post-metamorphic froglets. Faeces were collected
weekly and pooled by month; each metagenome is one monthly pool from one animal
(stems named `<animal>.<timepoint>`).

| cohort | metagenomes | animals | design |
|---|---|---|---|
| WF22 | 44 | 15 | three monthly pools per animal (one animal with two) |
| WF23 | 9 | 9 | one pool per animal |
| WF24 | 35 | 35 | one pool per animal |

**WF22.** Three egg masses; two *Basidiobolus* strains and control; three monthly
timepoints. Six of 50 metagenomes were excluded before cataloging: one control
with *Basidiobolus* reads, one that failed quality review, and four from animals
that died during the experiment.

**WF23.** Egg mass 3; 5 control and 4 UHM520.7734.

**WF24.** 36 animals sequenced across seven treatment codes; treated animals
received at least two inoculations within the month. One control with elevated
*Basidiobolus* reads was excluded, leaving 35 (32 from egg mass 4, 3 from egg
mass 2). One further animal was dropped from the experiment and is excluded from
treatment tests only, which therefore use 34 animals (4 controls, 5 per strain).
Tank is confounded with egg mass.

Because the WF24 exclusion was made after clustering, its proteins were removed
from the cluster tables and the support rule was reapplied; original
representatives were kept. Excluded samples are listed in
`metadata/wf24_excluded.tsv`, and affected numbers in `results/n88_number_table.tsv`.

| File | Purpose |
|---|---|
| `scripts/step25a_filter_n88.py` | Supported set, count matrix and sample table at 88 metagenomes |
| `scripts/step25b_membership_n88.py` | Cluster tables at 88 metagenomes |
| `scripts/step25c_number_table.py` | Old and new values of all reported numbers |

Output: `results/supported_reps_95_n88.txt`, `results/count_matrix_primary_n88.tsv.gz`,
`metadata/samples_88_treatment.tsv`, `results/clusters_{95,90,50}_n88.tsv`,
`results/n88_number_table.tsv`

### Step 1. Gene prediction

Proteins were predicted from assembled contigs with Prodigal v2.6.3 in metagenomic mode.

**Result.** 40,550,595 proteins from 89 metagenomes; 40,257,962 after excluding
one WF24 sample.

| File | Purpose |
|---|---|
| `scripts/step1c_prodigal_array.sh` | Gene prediction |
| `scripts/step1d_rename_45.sh` | Headers prefixed with source metagenome |
| `scripts/step1e_protein_counts.sh` | Protein counts per metagenome |

Output: `results/step1_protein_counts_45.tsv`

### Step 2. Clustering

Proteins were clustered with MMseqs2 at 95%, 90% and 50% amino-acid identity.

**Result.** 6,171,602 clusters at 95% identity, 5,352,641 at 90% and 2,918,630
families at 50%.

| File | Purpose |
|---|---|
| `scripts/step2_cluster.sh` | Database build and clustering at three thresholds |
| `scripts/step25b_membership_n88.py` | Cluster tables at 88 metagenomes |

Output: `catalog/db/LsPS_AA_{95,90,50}_cluster`, `results/clusters_{95,90,50}_n88.tsv`

### Step 3. Overlap with the earlier WF22 catalog

A random subsample of WF23 and WF24 proteins was searched with MMseqs2 against the
earlier WF22-only catalog, using query-side coverage. A positive control of WF22
proteins was recovered at 98.25%.

**Result.** 45.89% of WF23 and WF24 proteins match the earlier catalog at gene
level and 66.49% at family level.

| File | Purpose |
|---|---|
| `scripts/step2c_diagnostic_epyc.sh` | Subsample search, gene level |
| `scripts/step2d_control_and_tiers_n88.sh` | Positive control and family level |
| `scripts/step2e_covmode2_n88.sh` | Query-side coverage search |

### Step 4. Cluster support and cohort sharing

Each cluster was assigned the cohorts its members come from. Clusters with at
least three members from at least two metagenomes form the supported catalog.

**Result.** 55.29% of 95% clusters are singletons. The supported catalog has
2,067,011 clusters. Clusters shared by all three cohorts rise from 6.72% of all
clusters to 20.06% of supported clusters, and to 26.00% at family level. The top
1% of families hold 45.77% of all proteins.

| File | Purpose |
|---|---|
| `scripts/step3_cluster_composition_n88.sh` | Cohort composition of each cluster |
| `scripts/step4_nonsingleton_n88.sh` | Support filter and cluster size distribution |

### Step 5. Cohorts compared at equal sampling effort

Cohorts differ in metagenome count, so they were compared at 9 metagenomes each,
averaged over 10 random draws.

**Result.** WF22-exclusive clusters fall from 15.27% to 9.52% at equal effort;
WF24-exclusive clusters are 25.28%, and 24.33% are shared by all three cohorts.
Richness at 9 metagenomes is 765,122 clusters for WF22, 1,166,662 for WF23 and
2,080,538 for WF24, tracking per-sample protein yield rather than cohort.

| File | Purpose |
|---|---|
| `scripts/step5_equaln_rarefaction_n88.sh` | Equal-effort comparison and per-cohort accumulation |

### Step 6. Sequencing batch

Sequencing run, centre, platform and read depth were tabulated for all 88
metagenomes to check cohort comparisons against batch.

**Result.** WF23 and WF24 share one run and yield similar protein numbers per
metagenome (422,857 and 431,504), so cohort alone does not change yield. WF24
spans two runs that differ 1.75-fold (431,504 against 756,692) on 1.24-fold more
reads, driven by assembly size; gene density per Mb is constant. WF22 was
sequenced on separate runs, so WF22 against WF23 and WF24 cannot be separated
from batch.

| File | Purpose |
|---|---|
| `scripts/step6b_batch_table_n88.sh` | Batch table and run contrasts |

Output: `results/batch_table_88.tsv`

### Step 7. Rarefaction

Clusters were accumulated against proteins sampled rather than metagenomes, with
WF24 split by sequencing run so that cohort and run are not conflated.

**Result.** At 3,600,000 proteins, WF22 and WF23 give 1,169,677 and 1,159,262
clusters (0.9% apart). WF24 gives 1,295,364 and 1,687,560 in its two runs (30.3%
apart), while WF23 and WF24 within one run differ by 11.7%. Sequencing run, not
cohort, drives richness. The pooled curve reaches 6,171,602 clusters at
40,257,962 proteins and still adds 436,425 clusters between 35 and 40 million
proteins, so the catalog is not saturated.

| File | Purpose |
|---|---|
| `scripts/step7_rarefaction_n88.sh` | Accumulation against proteins sampled |
| `scripts/step7b_wf24_byrun_n88.sh` | Richness at matched depth, WF24 by run |

Output: `results/rarefaction_95_n88.tsv`

### Step 8. Annotation

Representatives at 95% identity were annotated with eggNOG-mapper and with Pfam
using hmmsearch.

**Result.** eggNOG annotates 2,980,052 of 6,171,602 representatives (48.29%) and
Pfam 2,709,544 (43.90%).

| File | Purpose |
|---|---|
| `scripts/step8f_eggnog_v2.sh` | eggNOG annotation |
| `scripts/step12b_pfam_hmmsearch.sh` | Pfam annotation |

### Step 9. Open reading frame completeness

Prodigal flags whether each protein has both a start and a stop codon or is
truncated at a contig edge. Completeness was summarized overall and by cluster size.

**Result.** 46.89% of the 40,257,962 proteins are complete. Among 95%
representatives, completeness rises with cluster size: 11.58% for singletons,
32.98% for pairs, 61.16% for three to five, 81.99% for six to twenty and 86.11%
for twenty-one or more. Singletons are mostly contig-edge fragments, which
supports using the 2,067,011 supported clusters as the catalog.

| File | Purpose |
|---|---|
| `scripts/step9_orf_completeness_n88.sh` | Completeness by cluster size |

Output: `results/orf_completeness_n88.tsv`

### Step 10. The unannotated fraction (eggNOG)

eggNOG annotation of the 6,171,602 representatives was summarized by cluster size
and completeness, then carried to the 50% family level with support filters.

**Result.** 48.29% of representatives are annotated, rising from 36.83% for
singletons to 90.53% for clusters of twenty-one or more. Proteins spanning a
whole short contig annotate better than complete open reading frames (65.1%
against 40.4% among singletons), so partial calls are not low-quality sequence.
At the family level, 76.03% of all families have no eggNOG-annotated member,
falling to 55.40% of the 1,020,958 supported families (565,611 families); the
value changes by less than 0.1 percentage points with a stricter support filter.

| File | Purpose |
|---|---|
| `scripts/step10_annotation_summary_n88.sh` | Annotated fraction by cluster size and completeness |
| `scripts/step11_interaction_and_supported_n88.sh` | Family-level unannotated fraction under support filters |

Output: `results/annotation_summary_n88.tsv`, `results/family_dark_fraction_n88.tsv`

### Step 11. Taxonomy of the annotated fraction

Broad clade was taken from eggNOG assignments of annotated representatives and
compared across cluster sizes.

**Result.** Of 2,980,052 annotated representatives, 92.64% are bacterial, 6.93%
eukaryotic, 0.27% archaeal and 0.16% viral. The most common groups are
Alphaproteobacteria (16.30%), Bacteroidetes (13.95%), Actinobacteria (13.34%),
Gammaproteobacteria (11.44%) and Betaproteobacteria (10.31%); Metazoa are 3.05%
and Fungi 2.57%. Eukaryotic share falls from 7.63% in singletons to 3.61% in
clusters of twenty-one or more, so host and dietary sequence is present but
sparse. The unannotated fraction has no taxonomic assignment from this analysis.

| File | Purpose |
|---|---|
| `scripts/step13_taxonomy_n88.sh` | Clade and assignment level by cluster size |

Output: `results/taxonomy_summary_n88.tsv`

### Step 12. Combined annotation and protein classes

eggNOG and Pfam (hmmsearch, gathering thresholds) were combined: a protein is
known if either source annotates it.

**Result.** The union annotates 3,222,677 of 6,171,602 representatives (52.22%).
Of the 2,067,011 supported representatives, 585,347 (28.32%) have no annotation.
At the 50% family level, the supported unannotated fraction is 51.47%, against
55.40% with eggNOG alone. Supported representatives fall into three classes:
K, with a Pfam domain (1,343,020, of which 76,515 carry only DUF or UPF domains);
KWP, eggNOG only (138,644); and U, unannotated (585,347).

| File | Purpose |
|---|---|
| `scripts/step12d_union_replevel_n88.sh` | eggNOG and Pfam union per representative |
| `scripts/step12e_family_union_n88.py` | Union at the 50% family level |
| `scripts/step12f_rep_classes_n88.sh` | K, KWP and U labels |

Output: `results/union_rep_level_n88.tsv`, `results/family_union_dark_n88.tsv`,
`results/rep_classes_n88.tsv`

### Step 13. Genomic and environmental unknowns

The 585,347 unannotated supported representatives were searched against NCBI
ClusteredNR (release 20260128, 470,748,714 sequences) with DIAMOND 2.1.24
--very-sensitive. A hit requires E <= 1e-10 with at least 50% coverage of both
query and subject.

**Result.** 430,299 representatives (20.82% of the supported catalog) have no hit
in ClusteredNR (environmental unknowns). 155,048 (7.50%) match only sequences
annotated as hypothetical (genomic unknowns). Across looser and stricter criteria
the environmental unknown share ranges from 63.92% to 76.63% of unannotated
representatives.

| File | Purpose |
|---|---|
| `scripts/step13d_extract_u.py` | Unannotated representatives to fasta |
| `scripts/step13f_nrclust_search.sh` | DIAMOND search against ClusteredNR |
| `scripts/step13h_gu_eu_split_n88.sh` | Genomic and environmental split under several criteria |
| `scripts/step13i_fourway_n88.py` | Final four-class labels |

Output: `results/step13_nrclust20260128_hits.tsv`, `results/fourway_classes_n88.tsv`

### Step 14. Structure of the unknown fraction

Unannotated representatives were screened for spurious open reading frames with
AntiFam v6.0 and clustered at 50% identity and 80% coverage.

**Result.** AntiFam flags 1,464 of 585,347 representatives (0.2501%). The
585,347 representatives form 536,030 families; 8,464 have three or more members
and 10 have 100 or more.

| File | Purpose |
|---|---|
| `scripts/step14b_antifam_n88.sh` | AntiFam screen |
| `scripts/step14c_novel_families_n88.sh` | Family sizes at 50% identity |

Output: `results/antifam_hits_n88.tblout`, `results/unk_clusters_50_n88.tsv`,
`results/unk_family_sizes_n88.tsv`

### Step 15. Abundance

Reads were mapped with bwa-mem2 2.3 to the 2,069,453 supported representatives
(nucleotide coding sequences). The matrix used in all analyses is restricted to
the 2,067,011 clusters supported at 88 metagenomes.

**Result.** Mean mapped fraction is 0.9587 in WF22, 0.9634 in WF23 and 0.9667 in
WF24, so the catalog represents all three years equally.

| File | Purpose |
|---|---|
| `scripts/step16h_supported_cds.py` | Mapping reference |
| `scripts/step16j_index_reference.sh` | Reference index |
| `scripts/step16q_map_array.sh` | Read mapping and counting |
| `scripts/step16u_matrix_n88.py` | Count matrix at 88 metagenomes |

Output: `results/count_matrix_primary_n88.tsv.gz`, `results/mapping_qc_88.tsv`

### Step 16. Treatment, egg mass and month: per-gene tests

Gene abundance was tested per gene with limma-voom (TMM normalization), each cohort
separately. WF22 samples are repeated measures of 15 animals, so animal was
blocked with duplicateCorrelation.

**Result.** Within-animal correlation is 0.2864; without blocking, WF22 treatment
contrasts return 5,005 and 7,263 genes, and with blocking 3 and 155. In the same
blocked model egg mass returns 57,120 and 51,808 genes and month up to 141,580.
Within the balanced egg mass 3 subset, treatment returns zero genes. WF23 returns
zero; WF24 (34 animals) returns zero for five of six strains and 3 genes for one.
The unannotated share of gene abundance does not differ by treatment in any cohort
(all p > 0.1).

| File | Purpose |
|---|---|
| `scripts/step17c_wf22_meta.sh` | WF22 design table |
| `scripts/step17g_limma_wf22_n88.R` | WF22, blocked and unblocked |
| `scripts/step17h_limma_wf2324_n88.R` | WF23 |
| `scripts/step17i_limma_wf24_n34_n88.R` | WF24, 34 animals |
| `scripts/step17b_class_shift_n88.py` | Unannotated abundance share by treatment |

Output: `results/limma_*_n88.csv`, `metadata/wf22_design.tsv`

### Step 17. Properties of the unknown fraction

Length, completeness, prevalence, persistence across years and abundance were
compared across four classes: K (Pfam domain), KWP (eggNOG only), GU (genomic
unknown) and EU (environmental unknown).

**Result.**
- Length: mean coding length is 1,045 bp for K, 795 for KWP, 498 for GU and 367
  for EU. Complete reading frames are 79.5%, 69.8%, 76.0% and 63.7%, so unknown
  genes are short but mostly complete.
- Prevalence: 82.77% of unannotated genes occur in more than 20 of 88 samples,
  3,496 in all 88, and 0.65% in five or fewer.
- Persistence: 496,266 unannotated genes occur in all three years and carry
  91.15% of unannotated reads; 18,896 single-year genes carry 1.23%, and 260
  genes receive no reads.
- Abundance: unannotated genes are 28.32% of supported genes, 7.04% of raw mapped
  reads and 19.32% after length normalization; the raw share is 7.11% in WF22,
  6.28% in WF23 and 7.16% in WF24 (per sample 4.38 to 16.12%).

| File | Purpose |
|---|---|
| `scripts/step16y_complete_by_class_n88.py` | Length and completeness by class |
| `scripts/step16v_prevalence_n88.py` | Prevalence by class |
| `scripts/step19a_core_unknown_n88.py` | Persistence across years |
| `scripts/step16w_dark_abundance_n88.py` | Abundance share, raw and length-normalized |
| `scripts/step18b_fig12.py` | Prevalence and abundance figures |

Output: `results/completeness_by_class_n88.tsv`, `results/prevalence_primary_n88.tsv.gz`,
`results/core_unknown_by_years_n88.tsv.gz`, `results/dark_abundance_by_sample_n88.tsv`

### Step 18. Taxonomic context of the unknown fraction

Each gene inherited the taxonomy of its contig, from existing UniRef50-based
contig assignments, to ask whether unannotated genes sit on less classifiable DNA.

**Result.** Contig taxonomy is available for 48 of 88 metagenomes, covering
888,183 of 2,067,011 genes (42.97%). Contigs are unclassified for 74.6% of K and
KWP genes, 74.4% of GU and 91.6% of EU genes. EU genes are 4.5% eukaryotic and
3.9% bacterial, against about 23% bacterial and 2% eukaryotic in the other
classes. Genes with no database match sit on DNA that is itself largely
unclassified, consistent with lineages lacking sequenced relatives.

| File | Purpose |
|---|---|
| `scripts/step20c_gene_tax_n88.py` | Contig taxonomy by gene and class |

Output: `results/gene_taxonomy_n88.tsv`

Correction (2026-10-16): the domain parse used here read only 'd_' lineage
tokens and missed '-_' labels. With the corrected parse
(scripts/step34c2c_taxv2.py, results/gene_taxonomy_n88_v2.tsv), EU genes on
classified contigs are mostly eukaryotic. Of genes with contig taxonomy
(results/step34c2c_summary.tsv), contigs are unclassified for 38.5% of K
(previously 74.6%), 25.9% of KWP (74.6%), 23.6% of GU (74.4%) and 33.8% of EU
genes (91.6%). EU genes are 53.8% eukaryotic (4.5%) and 12.3% bacterial
(3.9%); K genes are 56.0% bacterial and 5.5% eukaryotic, KWP 57.0% and 16.9%,
GU 67.5% and 8.1%. Among genes whose contig has a domain, EU is 81.2%
eukaryotic (93,809 of 115,463), against 8.9% for K, 22.8% for KWP and 10.7%
for GU. The statement above that genes with no database match sit on largely
unclassified DNA does not hold under the corrected parse.

### Step 19. Gene neighbourhood of unknown families

For the 8,464 unknown families with three or more members, annotated domains on
adjacent genes were counted. A domain beside at least half of a family's members
is taken as a conserved neighbour, a hypothesis about function rather than an
assignment.

**Result.** 3,614 families (42.7%) have a conserved annotated neighbour. The
strongest cases sit beside oxidoreduction and respiration domains (GSDH, Rieske,
COX1, adh_short, Aldo_ket_red), transporters (MFS_1, ATP_bind_1), regulators
(HTH_1, TetR-like, PhoU), surface proteins (Fimbrial, AsmA) and peptidases.
Members of these families are 6.29% of supported unannotated genes.

| File | Purpose |
|---|---|
| `scripts/step21b_neighbours_n88.py` | Neighbour domains of unknown families |

Output: `results/unk_family_neighbours_n88.tsv`

### Step 20. Community-level tests

Whole-community gene abundance (Bray-Curtis on Hellinger-transformed counts) and
gene presence (Jaccard, at least five reads) were tested with PERMANOVA. In WF22,
between-animal terms (treatment, egg mass) were tested by permuting whole animals
(42 samples, 14 animals) and on animal centroids (15 animals); month was permuted
within animals. A Pfam-defined carbohydrate subset (59,297 representatives, 7,891
tested) was analysed the same way and per gene.

**Result.**
- Treatment: R2 0.051 (p = 0.928) by whole-animal permutation and 0.069
  (p = 0.657) on centroids; presence R2 0.062 (p = 0.665) and 0.106 (p = 0.465).
  In WF24 (34 animals) R2 is 0.173 (p = 0.679) and 0.170 for presence
  (p = 0.770), at chance level for 6 of 33 degrees of freedom. Genes detected per
  sample do not differ (p = 0.3676).
- Egg mass: R2 0.327 by whole-animal permutation and 0.513 on centroids
  (both p = 0.001), with equal dispersions (p = 0.214). Month: R2 0.127 (p = 0.001).
- Carbohydrate subset: zero genes respond to treatment, against 1,541 and 2,046 for
  egg mass and 1,761 and 2,586 for month; composition gives treatment R2 0.053
  (p = 0.454) and egg mass 0.673 (p = 0.001).
- Treatment groups differ in dispersion (p = 0.003) but not in centroid.

| File | Purpose |
|---|---|
| `scripts/step22b_permanova_fix_n88.R` | Abundance composition, WF22 |
| `scripts/step22c_jaccard_n88.R` | Presence composition and richness, WF22 |
| `scripts/step22e_permanova_wf24_n34_n88.R` | Abundance and presence composition, WF24 |
| `scripts/step22d_cazy_subset_n88.R` | Carbohydrate subset |

Output: logs `step22b_n88` to `step22e_n88`, `results/limma_cazy_*_n88.csv`

### Step 21. Unknown proteins against NMPFams2

The unannotated supported representatives were filtered to match the NMPFams2
family set (Aplakidou et al. 2026): complete open reading frame and at least
35 aa, leaving 384,260 of 585,347 (65.65%). These were searched with MMseqs2
(`search -e 1e-3 --min-seq-id 0.3 -c 0.7 --cov-mode 0`, default -s 5.7) against
the representatives of NMPFams2 metagenome families (6,588,876) and isolate
families (2,229,233), each with at least 25 members. Coverage mode 0 (both
sequences) is used here rather than mode 2 because NMPFams2 families are built at
80% coverage of both sequences. A protein is represented if it has any hit.

**Result.** 1,787 of 266,848 environmental unknowns (0.67%) are represented in
either set (1,272 metagenome, 535 isolate, 20 both), so 99.33% match neither.
Genomic unknowns are represented far more often: 37,821 of 117,412 (32.21%),
mostly through isolate families (34,421). Over all 384,260 proteins, 39,608
(10.31%) are represented. Isolate families come from reference genomes and hold
annotated and unannotated proteins, so a match to one is not a match to an
uncharacterized family.

**Search sensitivity.** Searching representatives only underestimates matches. For
10,000 random queries (seed 26), 151 hit a metagenome family representative and
233 hit any member (1.51% against 2.33%); for isolate families 901 against 973.
In a positive control, 5,000 random non-representative NMPFams2 metagenome members
found their own family's representative in 4,873 cases (97.46%). Recovery stays
high for short queries (2,997 of 3,089 under 100 aa, 97.02%; control median
83 aa), and environmental unknowns are no shorter (median 94 aa), so their low
representation is not a short-query artifact.

| File | Purpose |
|---|---|
| `scripts/step26a_unk_linclust30.sh` | Complete-ORF and length filter |
| `scripts/step26b_get_nmpfams2.sh` | NMPFams2 download and md5 check |
| `scripts/step26b2_resume_iso.sh` | Resume isolate-family download |
| `scripts/step26b3_nmpfams2_search.sh` | Search against family representatives, 10,000-query sensitivity check |
| `scripts/step26b3_summary.py` | Representation of proteins and families |
| `scripts/step26d1_gu_eu_split.py` | Metagenome-family representation by genomic and environmental class |
| `scripts/step26d3_iso_gu_eu_split.py` | Isolate-family representation by class |
| `scripts/step26d4_union_gu_eu.py` | Union of metagenome and isolate representation |
| `scripts/step26e2_posctrl.sh` | Positive control |
| `scripts/step26e3_posctrl_by_length.py` | Positive control and query lengths by length bin |

Output: `results/unk_filtered_aplakidou.faa`, `results/nmpfams2_summary.tsv`,
`results/nmpfams2_{,iso_,union_}gu_eu.tsv`, `results/nmpfams2_posctrl_summary.tsv`,
`results/nmpfams2_posctrl_by_length.tsv`

### Step 22. Environmental-only 30% families and their recurrence

All 2,441,674 complete members of at least 35 aa in the 95% clusters of the
filtered representatives were clustered directly with MMseqs2 (`linclust
--min-seq-id 0.3 -c 0.8 --cov-mode 0`, as for NMPFams2). Each family was labelled
by the classes of its members' 95% representatives, and recurrence was counted
over metagenomes, animals and cohorts.

**Result.** The members form 412,858 families: 13,680 have at least 25 protein
records and 236 have at least 100, but counted as exact-unique sequences only 29
reach 25 and 1 reaches 100. 303,043 families are environmental-only, 109,260
genomic-only and 555 mixed. Of the environmental-only families, 740 contain at
least 3 distinct 95% representatives and 8 contain at least 10.

**Recurrence.** 30,780 environmental-only families occur in at least 10
metagenomes, but 7,338 of these reach 10 only through repeated sampling of WF22
animals, so recurrence is counted in animals. 69,605 families occur in at least
2 cohorts and 19,563 in all three. 23,442 occur in at least 10 animals, and 15,580
in at least 10 animals and at least 2 cohorts; this set is carried into Steps 24
and 25.

| File | Purpose |
|---|---|
| `scripts/step26a2_raw_counts.py` | Family sizes through 95% representatives, superseded by step26e1 |
| `scripts/step26d2_unique_counts.py`, `.sh` | Unique-sequence sizes through 95% representatives, superseded by step26e1 |
| `scripts/step26e1_direct30.sh` | Direct 30% clustering of member proteins |
| `scripts/step26f_direct30_gu_eu.py` | Genomic and environmental composition of families |
| `scripts/step26f2_eu_recurrence.py` | Metagenomes, animals and cohorts per environmental-only family |

Output: `results/unk_direct30_summary.tsv`, `results/unk_direct30_eu_families.tsv`,
`results/unk_direct30_eu_recurrence.tsv`, log `step26f`

### Step 23. *Basidiobolus* proteins in the catalog

Proteins from five *Basidiobolus* genomes (STP1710, STP1717, UHM207, UHM260,
UHM520) were searched against the catalog's 95% representatives with MMseqs2
(`search -e 1e-5 --min-seq-id 0.9 -c 0.8 --cov-mode 1`, coverage of the catalog
protein). The catalog protein in each best hit was called bacterial when its eggNOG
lineage or its contig's taxonomy was Bacteria. Full-length candidates were checked
for gene neighbours on their contig and for presence in control animals.

**Result.** 42,456 *Basidiobolus* proteins hit 14,943 catalog representatives. In
3,062 cases the best hit is a bacterial protein, but in 2,245 of these less than
half of the fungal protein is covered (a domain inside a longer protein). Requiring
80% coverage of both proteins leaves 338 queries and 61 catalog proteins. None of
the 61 occurs in a retained control sample; 3 occur in the excluded control
UHM586.41010, each also in 5 to 18 inoculated samples. 52 of the 61 are the only
gene on their contig and 58 have no annotated neighbour, so their origin is
unresolved. The one candidate with a bacterial neighbour lies on a 533 bp contig
that aligns end to end (blastn, E <= 1e-5) to all six *Basidiobolus* genomes,
including UHM516, which has no protein file and was not in the protein search
(27 copies in STP1710 at up to 94.65% identity), so the contig and its
"bacterial" neighbour are a multi-copy *Basidiobolus* repeat. There is no
convincing recent bacterial donor at >= 90% identity. Presence is from assembled proteins at 89 metagenomes, before the
WF24 exclusion.

| File | Purpose |
|---|---|
| `scripts/step27a_basid_vs_catalog.sh` | *Basidiobolus* proteins against catalog representatives |
| `scripts/step27b_label_hits.py` | Bacterial labels and best hits |
| `scripts/step27c_contig_context.py` | Contig neighbours of full-length candidates |
| `scripts/step27d1_candidate_presence.py` | Presence of candidates in control and inoculated samples |
| `scripts/step27d2_contig_blastn.sh` | blastn of the bacterial-context contig against the genomes |

Output: `results/basid_vs_catalog_labelled.tsv`, `results/basid_bacterial_candidates.tsv`,
`results/basid_candidates_context.tsv`, `results/basid_candidates_presence.tsv`,
`results/basid_contig_k141_15734_blastn.tsv`

### Step 24. Structural search of recurrent environmental-only families

Representatives of the 15,580 recurrent families (Step 22) were translated to 3Di
structural sequence with ProstT5 and searched with Foldseek (`search -e 1e-1`)
against PDB and CATH50. A hit is E <= 1e-3, set on a pilot of the 740 families with
at least 3 distinct 95% representatives, 740 length-matched annotated proteins and
the 740 pilot sequences reversed as decoys. ProstT5 databases carry no coordinates,
so hits are E-value based, with no TM-score.

**Calibration.** At E <= 1e-3 the pilot gives 7.57% hits in PDB and 7.84% in
CATH50, the length-matched annotated proteins 62.43% and 60.68%, and the reversed
decoys 0% in both.

**Result.** 218 of the 15,580 families (1.40%) have a hit in either database
(PDB 117, CATH50 196). The rate rises with length but stays below the pilot's in
every bin:

| length (aa) | recurrent set | pilot | annotated, length-matched | reversed |
|---|---|---|---|---|
| < 100 | 75/8,615 (0.87%) | 1/353 (0.28%) | 158/350 (45.14%) | 0/353 |
| 100-199 | 89/5,744 (1.55%) | 18/230 (7.83%) | 189/233 (81.12%) | 0/230 |
| 200-399 | 41/1,176 (3.49%) | 27/123 (21.95%) | 111/123 (90.24%) | 0/123 |
| >= 400 | 13/45 (28.89%) | 23/34 (67.65%) | 33/34 (97.06%) | 0/34 |
| all | 218/15,580 (1.40%) | 69/740 (9.32%) | 491/740 (66.35%) | 0/740 |

118 pilot families are also in the recurrent set and give the same 3 hits (2.54%)
in both runs, so the lower rate reflects family selection, not length or a
technical difference. The >= 400 aa bins are small (45 and 34), and no decoys were
run on the recurrent set.

| File | Purpose |
|---|---|
| `scripts/step28a_install_foldseek.sh` | Foldseek GPU build and ProstT5 weights |
| `scripts/step28a2_foldseek_epyc.sh` | ProstT5 weights download, rerun on epyc |
| `scripts/step28b_foldseek_dbs.sh` | PDB and CATH50 databases |
| `scripts/step28c_prostt5_gputest.sh` | ProstT5 GPU test |
| `scripts/step28d_inputs.py` | Pilot query set |
| `scripts/step28e_pilot_cpu.sh` | Pilot search, ProstT5 on CPU |
| `scripts/step28e2_convert.sh` | Pilot hit tables |
| `scripts/step28f_sets.py` | Length-matched annotated and reversed calibration sets |
| `scripts/step28f_calib.sh` | Calibration searches |
| `scripts/step28f_summary.py` | Calibration hit rates |
| `scripts/step28g_inputs.py` | Recurrent query set |
| `scripts/step28g_main.sh` | Recurrent-set search |
| `scripts/step28g_summary.py` | Recurrent-set hit rates |
| `scripts/step28g2_by_length.py` | Hit rates by length bin, pilot and recurrent overlap |

Output: `results/prostt5_calib_summary.tsv`, `results/prostt5_main_summary.tsv`,
`results/prostt5_by_length.tsv`

### Step 25. Predicted structures of environmental-only unknown families

ESMFold models were built for family representatives in three groups drawn from
the 15,580 environmental-only unknown families: 218 with a structural hit
(best E <= 1e-3 against PDB or CATH50), 200 with no hit at E <= 1e-1 in either
database, ranked by number of animals, then metagenomes, and 200 more no-hit
families of at least 100 aa not in the first set, ranked the same way. All 618
were folded with esm 2.0.0 (esmfold_3B_v1) on CPU only with `--chunk-size 64`;
618 of 618 completed, exit 0. A model counts as folded at mean pLDDT >= 70 and
pTM >= 0.5 unless it is helix-only (helix fraction >= 0.8, strand fraction < 0.05);
below either cutoff it is low confidence.

**Confidence by group.** HIT families (median 129.5 aa) have median pLDDT 62.8
and pTM 0.50, and 59 of 218 (27.1%) are folded. NOHIT (median 64.5 aa) has median
pLDDT 49.8 and pTM 0.25, with 5 of 200 (2.5%) folded; NOHIT_LONG (median 133.5 aa)
has 43.0 and 0.26, with 5 of 200 (2.5%) folded. Restricted to 100 aa and longer,
HIT gives 43 of 143 folded (30.1%, pLDDT 63.9, pTM 0.53), NOHIT 0 of 42 and
NOHIT_LONG 5 of 200, so the gap remains when short models are excluded. Most no-hit models are
low confidence (183 of 200 and 194 of 200).

**Folded no-hit models.** 10 no-hit models are folded (pLDDT 70.1 to 87.5, pTM
0.502 to 0.844). Four are compact, with radius of gyration 1.00 to 1.18 times the
value expected for a globular chain of that length (2.2 N^0.38 Å); the other six
are 1.36 or more. The four were checked on the Foldseek web server, which gave
the working labels DivIVA-like (99 aa), Cse2/CasB-like (161 aa), HTH-like (weak,
67 aa) and helical repeat-like (weak, 109 aa). No output was kept from those
searches; a local Foldseek search with TM-scores is pending. These are structural
resemblances, not novel folds.

**Presence.** Gene abundance (at least 5 reads on each family's supported 95%
representative, 88 metagenomes):

| family | metagenomes | animals | WF22 | WF23 | WF24 |
|---|---|---|---|---|---|
| DivIVA-like | 59 | 54 | 18/44 | 6/9 | 35/35 |
| Cse2/CasB-like | 50 | 35 | 28/44 | 6/9 | 16/35 |
| HTH-like (weak) | 75 | 55 | 31/44 | 9/9 | 35/35 |
| helical repeat-like (weak) | 40 | 24 | 30/44 | 5/9 | 5/35 |

WF22 was sequenced on separate runs (Step 6), so differences between WF22 and
the other cohorts cannot be separated from sequencing run.

**Gene neighbourhood.** Pfam domains on genes within three positions of every
family member were counted. The Cse2/CasB-like family sits in a type I-E cas
operon: Cse1 (CRISPR_Cse1) and Cas7 (Cas_CT1975) are within one gene in 37 of 37
members, and Cas5 (Cas_Cas5d) is within three genes in 35 of 37. The DivIVA-like
family's neighbours are UbiA (47 of 47) and DUF4229 (46 of 47), neither linked
to cell division, so that label is tentative. The HTH-like family sits beside
HTH_26 (73 of 73) and Uma2 (72 of 73).

| File | Purpose |
|---|---|
| `scripts/step29a_esmfold_input.py` | HIT, NOHIT and NOHIT_LONG representatives |
| `scripts/step29b_esmfold.sh` | ESMFold, CPU |
| `scripts/step29_peek_plddt.py` | Interim pLDDT check during the run |
| `scripts/step29_peek_v2.py` | Confidence, helix and strand fractions, compactness, folded call |
| `scripts/step29c_fold4_presence.py` | Assembly and read presence of the four compact families |
| `scripts/step29d_fold4_neighbours.py`, `.sh` | Neighbour domains of the four compact families |

Output: `results/esmfold_input_labels.tsv`, `results/esmfold_pdb/`,
`results/esmfold_peek_v2.tsv`, `results/step29_fold4_presence.tsv`,
`results/step29d_fold4_neighbours.tsv`, log `step29b_esmfold.29932169`

### Step 26. Local structural comparison of the folded no-hit models

The 10 folded no-hit ESMFold models from Step 25 were searched with a local Foldseek
(commit 463739e, run on CPU) against PDB and CATH50 with `--exhaustive-search 1
-e 10`, which reports TM-scores. For the four compact models, TM-scores were also
computed on Foldseek's own alignments from the saved web-server results (3Di/AA
mode), which include AlphaFold DB, Swiss-Prot and MGnify databases.

**Cse2/CasB-like (161 aa).** The best local hit in both databases is the CasB chain
of an *E. coli* Cascade structure (PDB 5h9e chain C; CATH50 5h9eC00): probability
1.000, qTM 0.618, tTM 0.654, query coverage 0.90, E 0.026 (PDB) and 0.020 (CATH50).
On the web results it aligns to Swiss-Prot Cse2 Q53VY0 at qTM 0.637, to AlphaFold DB
Cse2/CasB proteins at up to qTM 0.735, and to an uncharacterized MGnify protein at
qTM 0.833.

**The other three compact models.** None has a local hit with qTM >= 0.5 and
probability >= 0.5 except the HTH-like model, whose best CATH50 match (domain
1.20.58.110 from Q2FXY6) has qTM 0.531 at probability 0.569 and E 1.03. The
DivIVA-like model's best local match is CATH50 2wujA00 (qTM 0.359); on the web
results it aligns to AlphaFold DB DivIVA-domain proteins at qTM 0.798 but tTM
0.086, and to Wag31 at qTM 0.681, tTM 0.368, so it matches a DivIVA region only.
The helical repeat-like model's best match in either local database is 4jw3
(qTM 0.319). These remain structural resemblances, not novel folds. Three of the
six elongated models have matches with qTM >= 0.5 and probability >= 0.5, all at
E > 1.

| File | Purpose |
|---|---|
| `scripts/step29e_foldseek_local.sh` | Local Foldseek with TM-scores, PDB and CATH50 |
| `scripts/step29f_foldseek_web_tm.py` | TM-scores from the saved web-server results |

Output: `results/step29e_foldseek_summary.tsv`, `results/step29e_foldseek_{PDB,CATH50}.m8`,
`results/step29f_foldseek_web_tm.tsv`, `results/foldseek_web4/`

### Step 27. Gene neighbourhoods of the recurrent environmental-only families

Neighbourhoods of the 15,580 recurrent environmental-only families (Step 22) were
compared with 2,000 annotated families matched by the number of sequence variants.
A variant is a distinct member 95% representative; each variant contributes its
most common set of Pfam and eggNOG annotations within three genes. An annotation is
conserved when it occurs in at least 2 variants and at least 50% of a family's
variants, is at least 5 times more common than in the background of the same set,
and occurs near members of at least 5 families. Annotated families are 50% clusters
of Pfam-annotated proteins (seed 30). A first version that counted conservation per
animal (step30a) was replaced by this variant-based count because one strain
recurring across animals was counted as conservation.

**Most recurrent families are a single sequence variant.** 14,918 of 15,580
(95.75%) have one variant, so their recurrence is one sequence recovered in many
animals. The baseline was drawn to the same variant distribution (1,915 of its
2,000 families are single-variant).

**Multi-variant families are no more context-conserved than annotated ones.** Of
662 unknown families with at least 2 variants, 298 have an annotation conserved at
50% (285 at 75%, 281 at 90%), against 37 of 85 baseline families (30 and 30). The
share differs in both directions by variant bin:

| variants | unknown, conserved at 50% | baseline, conserved at 50% |
|---|---|---|
| 2 | 251/544 (46.1%) | 28/70 (40.0%) |
| 3 | 25/62 (40.3%) | 5/8 (62.5%) |
| 4 | 8/20 (40.0%) | 1/3 (33.3%) |
| 5-9 | 12/33 (36.4%) | 3/4 (75.0%) |
| >= 10 | 2/3 | 0/0 |

The baseline uses 50% clusters and the unknowns 30% families, so the comparison is
approximate.

**Eukaryotic neighbourhoods.** A family is flagged when Eukaryota is at least half
of its neighbour representatives with a classified contig domain. With the
corrected domain parse (Step 38), 12,322 of the 14,767 families with a classified
neighbour are flagged (`results/step34c2c_summary.tsv`); the parse used in this
step gave 3,030 of 4,631.

**Single-variant context.** Among the 14,918 single-variant families, the most
common neighbour set includes transcription or regulation annotations in 819,
metabolism 799, transport or membrane 714, phage or mobile elements 344, DNA
replication or repair 196, cell division or envelope 186, CRISPR-Cas or other
defence systems 84 and secretion 38. These are neighbour annotations of one
genome per family, not conserved calls. The Cse2/CasB-like family is one variant
(37 members, 26 animals, genus *Microbacterium*), whose most common neighbour set
holds CRISPR_Cse1 (1.00), Cas_CT1975 (1.00) and Cas_Cas5d (0.95): a canonical type
I-E neighbourhood of one sequence variant, not conservation across variants. The
15,580 families were later replaced by the strict_final families (Step 38),
because many were of fly origin (Steps 34-35).

| File | Purpose |
|---|---|
| `scripts/step30a_eu15580_neighbours.py`, `.sh` | Per-animal neighbourhood count, superseded |
| `scripts/step30b_eu15580_neighbours.py`, `.sh` | Variant-based neighbourhood conservation, matched baseline |

Output: `results/step30b_family_summary.tsv`, `results/step30b_baseline_comparison.tsv`,
`results/step30b_single_variant_context.tsv`, `results/step30b_family_annotations.tsv`,
`results/step30b_category_map.tsv`

### Step 28. Structures of defence-context single-variant families

The 84 single-variant families whose most common neighbour set holds a CRISPR-Cas,
restriction-modification or toxin-antitoxin annotation (Step 27) were folded with
ESMFold (11 models reused from Step 25, 73 new, same settings, terminal stop
character removed) and searched with the local Foldseek of Step 26. A match needs
qTM >= 0.5 and probability >= 0.5 in PDB or CATH50.

**Result.** 4 models are folded, 1 helix-only and 79 low confidence (median 69.5
aa, median pLDDT 49.45). Five families have a match, and only the Cse2/CasB-like
control matches the system of its neighbours (PDB 5h9e chain C, qTM 0.618, pLDDT
70.1, pTM 0.740). UHM581.41005__k141_82078_63, next to PemK_toxin, matches an
accessory protein of a relaxosome structure (PDB 5vfy chain D, E 4.5e-7, qTM 0.703,
tTM 0.824) on a low-confidence model (pLDDT 55.9, pTM 0.557).
UHM22.35746__k141_11288_106 (folded, next to Abi_2) matches 1r7j (qTM 0.545,
E 0.021); the two others match at E > 1. None of these four matches its neighbours.

| File | Purpose |
|---|---|
| `scripts/step31a_def84_input.py` | Defence-context family selection |
| `scripts/step31b_esmfold.sh` | ESMFold, CPU |
| `scripts/step31c_def84_summary.py`, `.sh` | Foldseek, confidence classes and summary |

Output: `results/step31_def84_summary.tsv`

### Step 29. Contig classification with Tiara

Contigs of at least 1 kb from all 88 assemblies were classified with Tiara 1.0.3
(`-m 1000 --probabilities`). Prokaryotic means bacteria, archaea or prokarya;
eukaryotic means eukarya or organelle. The excluded control UHM586.41010 was run
separately and used only to label representatives whose own contig is in that
assembly. Reads were not depleted of host DNA before assembly.

**Result.** 3,795,590 contigs were classified: 3,147,020 prokaryotic (831,604 of
them >= 3 kb), 630,389 eukaryotic (196,589 >= 3 kb) and 18,181 unknown.

| File | Purpose |
|---|---|
| `scripts/step32a_tiara_array.sh` | Tiara on the 88 assemblies |
| `scripts/step32a2_tiara_41010.sh` | Tiara on UHM586.41010 |

Output: `results/tiara/<stem>.tiara.tsv`

### Step 30. Tiara labels of protein classes and families

Each supported representative takes the Tiara class of its own contig, in bands of
>= 3 kb, 1-3 kb and < 1 kb (not classified).

**Classes.** Within the >= 3 kb band, environmental-only unknowns (EU) are mostly on
eukaryotic contigs, the other classes mostly on prokaryotic ones:

| class | prokaryotic | eukaryotic | unknown | share on contigs < 1 kb |
|---|---|---|---|---|
| K | 95.36% | 4.54% | 0.10% | 11.83% |
| KWP | 89.02% | 10.78% | 0.20% | 22.61% |
| GU | 92.49% | 7.28% | 0.23% | 25.52% |
| EU | 14.09% | 84.95% | 0.96% | 33.22% |

In the 1-3 kb band EU is 19.75% prokaryotic and 78.90% eukaryotic.

**Families.** Of the 15,580 recurrent families (233,120 members), 12,542 (80.50%)
have most members on eukaryotic contigs >= 3 kb and 2,895 (18.58%) on prokaryotic
ones. The Cse2/CasB-like family has all 37 members on prokaryotic contigs.

**NMPFams2.** EU representatives on prokaryotic contigs are also mostly outside
NMPFams2: 23,781 of 24,415 (97.40%) unmatched at >= 3 kb and 8,341 of 8,532
(97.76%) at 1-3 kb, against 201,646 of 202,198 (99.73%) on eukaryotic contigs and
265,061 of 266,848 (99.33%) overall. The eukaryotic source is identified in Steps
31-35: mostly the feeder fly, *Drosophila*, with little frog DNA.

| File | Purpose |
|---|---|
| `scripts/step32b_tiara_map.py`, `.sh` | Tiara class per representative, family and NMPFams2 group |

Output: `results/step32b_class_tiara.tsv`, `results/step32b_family_tiara.tsv`,
`results/step32b_nmpfams2_tiara.tsv`

### Step 31. Wood frog DNA in the assemblies

Contigs were mapped to the wood frog genome (*Lithobates sylvaticus*,
GCA_028564925.1, chromosome level; `catalog/host/SOURCE.txt`) with minimap2 2.30
(`-x asm10 --secondary=no -c`). A contig is host when its merged aligned length
covers at least 50% of it at at least 95% identity (matches over alignment block
length from the CIGAR alignment); 90% is the sensitivity threshold. Four contig sets
were tested: A, Tiara-eukaryotic contigs >= 1 kb; B, contigs carrying an EU
supported representative; F, contigs carrying a member of the 15,580 recurrent
families; C, a 2% random sample of Tiara-prokaryotic contigs >= 3 kb (seed 33) as
a negative control.

**Positive control.** 5,000 fragments simulated from the frog genome (1%
substitutions, lengths drawn from set A, seed 33) are called host in 4,984 cases
(99.68%) at >= 95% and 4,990 (99.80%) at >= 90%; median identity 0.9902.

**Result.** Frog DNA is a minor component:

| set | contigs | host at >= 95% | host at >= 90% |
|---|---|---|---|
| A, Tiara-eukaryotic >= 1 kb | 630,389 | 3,235 (0.51%) | 3,474 (0.55%) |
| B, contigs with an EU representative | 323,147 | 1,327 (0.41%) | 1,405 (0.43%) |
| F, contigs with a recurrent-family member | 134,321 | 194 (0.14%) | 213 (0.16%) |
| C, negative control | 16,632 | 8 (0.05%) | 16 (0.10%) |

Of the 430,299 EU representatives, 1,341 (0.31%) are on host contigs, most of them
on contigs < 1 kb (1,084 of 142,930, 0.76%; at most 0.02% on contigs >= 3 kb).
14 of the 15,580 families have at least half of their members on host contigs; the
Cse2/CasB-like family has none. A first mapping run without `-c` (step33c v1)
reported only seed-chain matches as identity (median 0.777 on simulated frog
fragments) and was replaced by this run; its files are marked superseded.

| File | Purpose |
|---|---|
| `scripts/step33b_host_index.sh` | minimap2 index of the frog genome |
| `scripts/step33c_prep.py`, `.sh` | Contig sets A, B, C and F |
| `scripts/step33c_map_array.sh`, `scripts/step33c_summary.py`, `.sh` | First mapping run, superseded |
| `scripts/step33c2_map_array.sh` | Mapping with `-c` |
| `scripts/step33c2_summary.py`, `.sh` | Host calls per contig, EU representative and family |
| `scripts/step33d_sim.py`, `scripts/step33d_poscontrol.sh`, `scripts/step33d_summary.py` | Simulated positive control, first run |
| `scripts/step33d2_poscontrol.sh`, `scripts/step33d2_summary.py` | Positive control with `-c` |

Output: `results/step33c_contig_host_v2.tsv`, `results/step33c_eu_host_v2.tsv`,
`results/step33c_family_host_v2.tsv`, `results/step33d_sim_v2/`

### Step 32. Tiara against read and contig taxonomy

The Tiara-eukaryotic fraction of contigs >= 3 kb per sample was compared with
Kraken2 (PlusPFP) read classifications, and Tiara classes with the UniRef50 LCA
domain of the same contigs. Cohort is confounded with sequencing run, so this is
descriptive.

**Per sample.** The Tiara-eukaryotic fraction of contigs >= 3 kb has a median of
0.95% (0.13% to 62.9%); 27 samples are at 10% or more, and 10 samples hold 60.7% of
all Tiara-eukaryotic contigs >= 3 kb (196,589 of 1,031,343 pooled, 19.06%).
Spearman rho with the Kraken2 Eukaryota read percentage is 0.53 over 88 samples
(WF22 0.44, WF23 0.60, WF24 0.81). In PlusPFP the only chordate with reads is
*Homo sapiens*, so frog DNA is not visible to Kraken2. WF22 holds repeated samples
of 15 animals, so rho is given without a p-value.

**Contig domain.** With the corrected domain parse (Step 38), Tiara and the LCA
domain agree on 113,168 of 114,528 contigs >= 3 kb with both a Tiara call and a
domain (98.81%); 34,755 of 36,386 Eukaryota-labelled contigs (95.52%) are
Tiara-eukaryotic.

| File | Purpose |
|---|---|
| `scripts/step33e_tiara_vs_kraken.py`, `.sh` | Per-sample comparison and Tiara against LCA domain |

Output: `results/step33e_tiara_vs_kraken.tsv`, `results/step33e_tiara_vs_lca.tsv`;
corrected-parse values in `results/step34c2c_summary.tsv`

### Step 33. 18S rRNA genes on Tiara-eukaryotic contigs

18S rRNA genes were predicted with barrnap 0.9 (`--kingdom euk`) on
Tiara-eukaryotic contigs >= 1 kb, and hits >= 500 bp were assigned to their top
SILVA 138.1 NR99 match with vsearch 2.21.1 (`--id 0.70`). The counts are numbers of
18S-bearing contigs, not abundances.

**Result.** 729 18S genes were classified, at least one in each of the 88 samples.
Insecta account for 96 (in 86 samples), 92 of them *Drosophila melanogaster*
(median identity 99.6%). Fungi account for 272 (86 samples), mostly yeasts
(*Pichia occidentalis* 75, *Graphium penicillioides* 66, *[Candida]
stellimalicola* 27, *Galactomyces candidus* 18); *Basidiobolus* gives 37 in 34
samples (32 *B. ranarum*). Mites typical of fly cultures (*Sancassania* 18,
*Hormosianoetus* 12) and plants or algae (94) are also present. Amphibia give 10
contigs in 9 samples (all *Rana catesbeiana* top hits, 99.9%), and no Ciliophora
or Opalinata 18S was found. 39 top hits to *Solanum chacoense* are at 87.1% to
90.6% identity and are not potato.

| File | Purpose |
|---|---|
| `scripts/step33f_barrnap_array.sh` | barrnap on Tiara-eukaryotic contigs |
| `scripts/step33f_collect.py` | 18S genes >= 500 bp |
| `scripts/step33f_classify.sh` | vsearch against SILVA 138.1 |
| `scripts/step33f_summary.py` | Composition by group and cohort |

Output: `results/step33f_18S_classified.tsv`, `results/step33f_composition.tsv`

### Step 34. *Drosophila* DNA in the assemblies

The four contig sets of Step 31 were mapped to *D. melanogaster* dm6
(GCF_000001215.4; `catalog/fly/SOURCE.txt`) with the same settings and rule. WF22
froglets were fed wingless fruit flies (*Drosophila* sp.), according to the draft
Methods of the experiment (Froglet Husbandry); feeding in WF23 and WF24 is not
recorded in the project files.

**Positive control.** 5,000 simulated fly fragments are called fly in all 5,000
cases at >= 95%; median identity 0.9902.

**Result.** Fly DNA is the largest identified eukaryotic component:

| set | contigs | fly at >= 95% | fly at >= 90% |
|---|---|---|---|
| A, Tiara-eukaryotic >= 1 kb | 630,389 | 269,721 (42.79%) | 277,510 (44.02%) |
| B, contigs with an EU representative | 323,147 | 102,598 (31.75%) | 106,215 (32.87%) |
| F, contigs with a recurrent-family member | 134,321 | 79,966 (59.53%) | 82,115 (61.13%) |
| C, Tiara-prokaryotic sample | 16,632 | 70 (0.42%) | 70 (0.42%) |

141,359 of the 430,299 EU representatives (32.85%) are on fly contigs, and 11,962
of the 15,580 recurrent families (76.78%) have at least half of their members on
fly contigs (12,232, 78.51%, at >= 90%). Fly and frog calls overlap on 1 contig.
The Cse2/CasB-like family has no fly member.

| File | Purpose |
|---|---|
| `scripts/step33g_sim.py`, `scripts/step33g_control.sh`, `scripts/step33g_sim_summary.py` | Simulated positive control |
| `scripts/step33g_map_array.sh` | Mapping to dm6 |
| `scripts/step33g_summary.py`, `.sh` | Fly calls per contig and family |

Output: `results/step33g_contig_fly.tsv`, `results/step33g_family_fly.tsv`,
`results/step33g_sim/`

### Step 35. Placement of fly-derived EU representatives on fly genes

Each EU representative on a fly contig was projected through the contig's
alignment onto dm6 and classed by the RefSeq annotation it overlaps (at least
30 bp, first match wins), then searched against the dm6 proteome with DIAMOND
blastp `--sensitive` (E <= 1e-3).

**Result.** Most fly-derived EU representatives are open reading frames outside
fly coding sequence in that frame:

| placement | EU representatives on fly contigs | with a dm6 protein hit |
|---|---|---|
| CDS, same strand, in frame | 6,697 (4.74%) | 95.46% (median identity 98.6%) |
| CDS, same strand, other frame | 2,936 (2.08%) | 5.76% |
| CDS, opposite strand | 14,437 (10.21%) | 0.57% |
| UTR or non-coding exon | 14,659 (10.37%) | 0.40% |
| intron | 57,882 (40.95%) | 0.06% |
| intergenic | 43,586 (30.83%) | 0.18% |
| not projected | 1,162 (0.82%) | 2.32% |
| all | 141,359 | 4.84% |

For the representatives of the fly-majority families (10,878 on a fly contig),
11.85% are in frame with fly CDS and 38.21% in introns. At least half
soft-masked: 5.98% of all and 1.99% of the family representatives.

| File | Purpose |
|---|---|
| `scripts/step33h_classify.py` | Projection onto dm6 and annotation class |
| `scripts/step33h_summary.py` | DIAMOND, repeat and class summary |
| `scripts/step33h.sh` | Job script |

Output: `results/step33h_rep_class.tsv`, `results/step33h_diamond_dm6.tsv`,
`results/step33h_summary.tsv`

### Step 36. Origin label for every supported representative

Each of the 2,067,011 supported representatives was labelled from its own contig,
first match in this order: fly, frog (Steps 31 and 34; contigs not mapped before
were mapped with the same settings), Tiara-eukaryotic other, Tiara-prokaryotic,
Tiara-unknown, or < 1 kb unclassified. The strict set is Tiara-prokaryotic contigs
>= 1 kb that are not fly or frog; the broad set adds Tiara-unknown and < 1 kb
contigs and is used only as a sensitivity comparison.

**Result.** EU representatives are mostly of eukaryotic or unclassified origin:

| origin | K | EU | all |
|---|---|---|---|
| fly | 20,117 (1.50%) | 141,359 (32.85%) | 175,931 (8.51%) |
| frog | 1,117 (0.08%) | 1,340 (0.31%) | 3,295 (0.16%) |
| Tiara-eukaryotic other | 45,358 (3.38%) | 146,214 (33.98%) | 203,489 (9.84%) |
| Tiara-prokaryotic (strict) | 1,120,953 (83.47%) | 41,392 (9.62%) | 1,362,125 (65.90%) |
| Tiara-unknown | 1,415 (0.11%) | 1,735 (0.40%) | 3,674 (0.18%) |
| < 1 kb unclassified | 154,060 (11.47%) | 98,259 (22.84%) | 318,497 (15.41%) |

On contigs < 1 kb, 30.50% of EU representatives are on fly contigs. The strict set
holds 1,362,125 representatives (EU 41,392, 3.04%), the broad set 1,684,296 (EU
141,386). The strict set was refined further in Step 38.

| File | Purpose |
|---|---|
| `scripts/step34a_prep.py` | Representative contigs not mapped before |
| `scripts/step34a_map_array.sh` | Frog and fly mapping of those contigs |
| `scripts/step34a_origin.py`, `.sh` | Origin label per representative |

Output: `results/step34a_rep_origin.tsv`, `results/step34a_origin_by_class.tsv`

### Step 37. Catalog statistics of the final microbial set

Statistics were recomputed for the final set (strict_final, Step 38), with the
earlier sets beside them. Read shares are percentages of all reads on the
2,067,011 supported representatives, raw and normalised by CDS length; a gene is
present in a sample at >= 1 read, and at >= 5 reads as a stricter criterion.

**Result.** strict_final holds 1,359,273 representatives and 92.99% of mapped
reads (83.94% length-normalised):

| class | representatives | mean CDS (bp) | complete ORF | mean samples present (>= 1 / >= 5 reads) | in all three cohorts (>= 1 / >= 5 reads) | reads (raw / length-normalised) |
|---|---|---|---|---|---|---|
| K | 1,118,691 (82.30%) | 1,121 | 88.09% | 52.61 / 33.67 | 95.19% / 64.71% | 85.73% / 70.83% |
| KWP | 94,261 (6.93%) | 901 | 85.89% | 45.79 / 27.73 | 91.20% / 53.71% | 3.84% / 4.36% |
| GU | 105,207 (7.74%) | 523 | 91.40% | 37.95 / 22.64 | 80.81% / 39.39% | 2.64% / 5.83% |
| EU | 41,114 (3.02%) | 372 | 78.12% | 37.05 / 22.72 | 72.13% / 35.17% | 0.79% / 2.92% |

In the whole catalog EU was 430,299 representatives (20.82%) with 3.56% of reads
(11.07% length-normalised); in the broad set 141,386 (8.39%) with 1.87%. Presence in
all three cohorts falls sharply from >= 1 to >= 5 reads in every class, so it
depends on sequencing depth. Per sample, reads on fly-derived representatives are a
median 0.47% of mapped reads (0.06% to 9.09%; WF22 0.52%, WF23 0.53%, WF24 0.45%);
cohort is confounded with sequencing run.

| File | Purpose |
|---|---|
| `scripts/step34b_catalog_stats.py`, `.sh` | Strict, broad and whole-catalog statistics; fly read share |
| `scripts/step34b2_final_stats.py`, `.sh` | strict_final statistics; NMPFams2 query filter |

Output: `results/step34b2_catalog_stats_final.tsv`, `results/step34b_catalog_stats.tsv`,
`results/step34b_fly_read_share.tsv`

### Step 38. Environmental-only families of the final microbial set

The strict EU representatives (41,392) were clustered at 30% identity with
`mmseqs linclust --min-seq-id 0.3 -c 0.8 --cov-mode 0` (as in Step 22, but on the
representatives, terminal stop removed). A family is recurrent when its members
come from at least 10 animals and at least 2 cohorts. Three refinements followed.

**Genus check.** Majority genus of each family's member contigs from the UniRef50
LCA (48 of 88 samples have contig taxonomy): 2,090 of the 4,561 strict recurrent
families carry a genus label. Fly-associated bacterial genera are the majority
genus of 14 (0.67%; Providencia 13, Leuconostoc 1), against 0.48% of non-recurrent
families and 0.39% of 5,000 random strict K representatives. The commonest are
*Sphingobacterium* 184, *Microbacterium* 165, *Enterococcus* 84, *Flavobacterium*
77 and *Corynebacterium* 75. 22 families had *Drosophila* as majority genus; their
23 representative contigs failed the fly rule on identity (19), coverage (3) or
lacked an alignment (1), and Tiara had called all of them prokaryotic.

**Domain parse.** The contig domain had been read only from `d_` lineage tokens,
which missed the `-_` form used in most samples with contig taxonomy; the
corrected parse reproduces the old file row for row and moves 376,459 genes out of
"unclassified" (see the note at Step 18).

**strict_v2 and strict_final.** strict_v2 removes representatives on contigs that
are fly at >= 90% identity (same coverage rule) or carry a Eukaryota LCA lineage:
88, 2,616 and 143 by fly only, lineage only and both (EU 57, 150 and 66).
strict_final then excludes the 5 single-representative recurrent families still
dominated by *Drosophila* (each failed the fly mapping rule). strict_final holds
1,359,273 representatives: K 1,118,691, KWP 94,261, GU 105,207 and EU 41,114.

**Result on strict_final.** 39,926 families, 38,916 (97.47%) with one
representative. 4,509 are recurrent (5,107 representatives), and 4,041 of them
(89.62%) are a single sequence variant. Of the strict_final EU representatives,
30,070 are complete ORFs of >= 35 aa and were searched against NMPFams2 (Step 21;
8,996 partial, 2,710 shorter, 662 both); 29,267 of them (97.33%) are unmatched,
and 3,868 of 4,127 (93.72%) among recurrent-family representatives. The
Cse2/CasB-like representative (UHM45.35758__k141_100597_236) stays in a recurrent
family of one variant: 45 member proteins from 32 animals in all three cohorts.
These families replace the 15,580 recurrent families of Step 22, which were built
before eukaryotic and diet-derived contigs were removed.

| File | Purpose |
|---|---|
| `scripts/step34c_prep.py`, `scripts/step34c_strict_eu.sh` | Strict EU sequences, 30% clustering |
| `scripts/step34c_families.py` | Recurrence, NMPFams2, Cse2 control |
| `scripts/step34c2_genus.py`, `.sh` | Majority genus of families and baselines |
| `scripts/step34c2b_drosophila.py`, `.sh` | Drosophila-majority families, contig details |
| `scripts/step34c2c_taxv2.py` | Corrected domain parse and the summaries that use it |
| `scripts/step34c2c_strict_v2.py`, `scripts/step34c2c.sh` | strict_v2 and family statistics |
| `scripts/step34c2c_exclude.py`, `.sh` | Exclusion of 5 families; strict_final |

Output: `results/step34c_strict_eu_families.tsv`, `results/step34c_strict_eu_nmpfams2.tsv`,
`results/step34c2_recurrent_genus.tsv`, `results/step34c2b_drosophila_families.tsv`,
`results/gene_taxonomy_n88_v2.tsv`, `results/step34c2c_summary.tsv`,
`results/step34c2c_strict_v2_families.tsv`, `results/step34c2c_excluded_families.tsv`,
`results/step34a_rep_origin_v2_final.tsv`, `results/step34b2_nmpfams2_filter.tsv`

### Step 39. Egg-mass recurrence of environmental-only genes

Presence of each strict_final EU representative was scored per animal (any member
protein assembled from the animal's samples; >= 5 reads as a sensitivity), with
41,114 random strict_final K representatives as a baseline (seed 343). WF22 uses
months 2 and 3, the months all 15 animals share (74 of 88 samples retained). A
representative is egg-mass restricted when it is present in at least half and at
least 3 of an egg mass's animals and absent from all animals of the other egg
masses that year. The null permutes egg-mass labels among animals within the year
(1,000 permutations, group sizes kept). WF23 has one egg mass; in WF24 egg mass
equals tank and egg mass 2 has 3 animals.

**Result.** In WF22, 948 EU representatives are egg-mass restricted against a null
mean of 67.5 (p = 0.001), and 616 K representatives against 77.6 (p = 0.001). With
the same permutations applied to both sets, EU exceeds K by 332 (null mean -7.1,
sd 45.4, p = 0.001); at >= 5 reads the difference is 454 (p = 0.001). Most come
from egg mass 3 (EU 837, K 576). WF24 egg mass 4 is not above null (EU p = 0.10).

| WF22 egg mass | animals | EU restricted (null) | K restricted (null) |
|---|---|---|---|
| 1 | 3 | 48 (4.9) | 29 (4.2) |
| 2 | 3 | 63 (5.8) | 11 (3.9) |
| 3 | 9 | 837 (56.8) | 576 (69.6) |

**Egg mass, not inoculation.** STP1710.7 occurs only in egg mass 3. On average,
each STP1710.7 animal carries 68.9% of the egg-mass-3-restricted EU representatives
(control 74.6%, STP1717.1 84.8%), and none of them is confined to one treatment
group. Egg masses were reared in separate tanks
(Dallas et al. manuscript draft, 22UHMwf_draft_V1,
Methods: Tadpole Husbandry), so egg mass means clutch plus rearing tank.

| File | Purpose |
|---|---|
| `scripts/step34c3_eggmass.py`, `.sh` | Egg-mass recurrence and permutation null |
| `scripts/step34c3b_checks.py`, `.sh` | Paired EU against K test; treatment check within egg mass 3 |

Output: `results/step34c3_eggmass_recurrence.tsv`, `results/step34c3b_checks.tsv`

### Step 40. Structural search of the final recurrent families

One representative per strict_final recurrent family (4,509, median 80 aa) was
searched as in Step 24, with the terminal stop removed: ProstT5 3Di, Foldseek
`search -e 1e-1` against PDB and CATH50, hit at best E <= 1e-3. Controls are 4,509
length-matched strict_final K representatives (10-aa bins, seed 34; where a bin had
too few K representatives, the rest came from the next bin up) and the 4,509
queries reversed.

**Result.** 186 of 4,509 families (4.13%) have a hit, against 50.21% of the
length-matched K representatives and 0.33% of the reversed decoys, so about 15 of
the 186 are expected to be false:

| length (aa) | recurrent families | K, length-matched | reversed |
|---|---|---|---|
| < 100 | 9/2,722 (0.33%) | 791/2,722 (29.06%) | 0/2,722 |
| 100-199 | 41/1,154 (3.55%) | 878/1,154 (76.08%) | 5/1,154 (0.43%) |
| 200-399 | 91/509 (17.88%) | 474/509 (93.12%) | 2/509 (0.39%) |
| >= 400 | 45/124 (36.29%) | 121/124 (97.58%) | 8/124 (6.45%) |
| all | 186/4,509 (4.13%) | 2,264/4,509 (50.21%) | 15/4,509 (0.33%) |

56 of the families with a hit in Step 24 have the same sequence among these
queries, and all 56 hit again. The Cse2/CasB-like representative reaches E 0.041
(CATH50 4h7aB00), short of the threshold; its CasB match was found by ESMFold and
Foldseek with TM-scores (Steps 25 and 26).

| File | Purpose |
|---|---|
| `scripts/step34d_inputs.py` | Query, length-matched K and reversed sets |
| `scripts/step34d_prostt5.sh` | ProstT5 and Foldseek, CPU |
| `scripts/step34d_summary.py` | Hit rates by length, decoy rate, overlap with Step 24 |
| `scripts/step34d2_kfix_inputs.py`, `step34d2_kfix_merge.py`, `step34d2_kfix.sh` | Replacement of 225 duplicated K controls |

Output: `results/step34d2_by_length_fixed.tsv`, `results/step34d_query_hits.tsv`,
`results/step34d_query_map.tsv`

### Step 41. Treatment, egg mass and month on the final set

The per-gene (Step 16) and community-level (Step 20) tests were rerun unchanged on
strict_final genes and on strict_final EU genes alone, with TMM normalisation
computed within each gene set. Added: a WF23 community test (126 distinct
labelings), month within animal for presence, presence dispersion, and WF24 with
35 animals as a sensitivity (primary 34, as before). Egg mass and treatment are
constant within WF22 animals and are untestable by within-animal permutation.

**Community.** All strict_final genes (EU genes in brackets):

| term | test | R2 | p |
|---|---|---|---|
| egg mass, abundance | whole-animal permutation, 42 samples | 0.342 (0.321) | 0.001 (0.001) |
| egg mass, abundance | animal centroids, 15 | 0.516 (0.491) | 0.001 (0.001) |
| month, abundance | within animal, 44 samples | 0.118 (0.107) | 0.001 (0.001) |
| treatment, abundance | whole-animal permutation | 0.052 (0.054) | 0.926 (0.877) |
| treatment, abundance | animal centroids | 0.069 (0.075) | 0.665 (0.591) |
| egg mass, presence | whole-animal permutation | 0.187 (0.235) | 0.001 (0.001) |
| egg mass, presence | animal centroids | 0.371 (0.443) | 0.001 (0.001) |
| month, presence | within animal | 0.102 (0.083) | 0.001 (0.001) |
| treatment, presence | whole-animal permutation | 0.071 (0.071) | 0.617 (0.644) |
| WF23 treatment, abundance | free permutation, 9 | 0.147 (0.154) | 0.245 (0.211) |
| WF24 treatment, abundance | free permutation, 34 | 0.175 (0.182) | 0.614 (0.494) |

Treatment groups differ in dispersion (p = 0.005; EU 0.004) but not in centroid.

**Per gene.** limma-voom blocked on animal (within-animal correlation 0.3374; EU
0.3059), genes at adjusted p < 0.05 among 294,643 tested (EU 5,604): treatment
STP1710.7 2 (0) and STP1717.1 0 (10); egg mass 2 and 3 48,380 and 46,098 (761
and 822); month 2 and 3 90,352 and 111,320 (1,564 and 1,660). STP1710.7 within egg
mass 3 gives 0 (1). WF23 gives 0, and WF24 0 for all six strains with 34 or 35
animals. Without blocking, treatment gives 4,281 and 7,884 genes. Egg mass and
month structure gene abundance and presence; treatment does not, and the whole-
catalog results of Steps 16 and 20 hold on the final set.

| File | Purpose |
|---|---|
| `scripts/step34e_counts.py` | Count matrices of the two gene sets |
| `scripts/step34e_hostfactor.R` | Community and per-gene tests for one gene set |
| `scripts/step34e_summary.py` | New values beside the whole-catalog values |
| `scripts/step34e_hostfactor.sh` | Job script |

Output: `results/step34e_summary.tsv`, `results/step34e_strict_final_permanova.tsv`,
`results/step34e_strict_final_limma.tsv`, `results/step34e_strict_final_EU_permanova.tsv`,
`results/step34e_strict_final_EU_limma.tsv`, `results/step34e_limma/`

## Software

Prodigal V2.6.3, MMseqs2 13-45111, eggNOG-mapper 2.1.9 (eggNOG 5.0.2),
HMMER 3.4 (Pfam), BBTools read statistics, Python 3.9 with NumPy 1.26.4.

## Repository layout

    scripts/   analysis and submission scripts
    metadata/  sample lists and per-sample metadata
    results/   summary tables
    catalog/   protein catalog and databases, not tracked
    logs/      job logs, not tracked
