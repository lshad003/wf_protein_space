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

## Software

Prodigal V2.6.3, MMseqs2 13-45111, eggNOG-mapper 2.1.9 (eggNOG 5.0.2),
HMMER 3.4 (Pfam), BBTools read statistics, Python 3.9 with NumPy 1.26.4.

## Repository layout

    scripts/   analysis and submission scripts
    metadata/  sample lists and per-sample metadata
    results/   summary tables
    catalog/   protein catalog and databases, not tracked
    logs/      job logs, not tracked

## Conventions

No count is stated without the output file that produced it. Coverage mode is stated explicitly for every search, since
it changes recovery by more than thirty percentage points. Cluster-level claims
are made on support-filtered clusters, and cohort comparisons at matched
sampling effort.

## Interpretation rules

1. Catalog richness is interpreted only after matching protein depth or
   accounting for sequencing run.
2. Cohort-exclusive clusters are interpreted only after support filtering, since
   singletons are necessarily cohort-exclusive and are dominated by contig-edge
   fragments.
3. Egg mass is tested only within cohorts where the sequenced samples provide
   replication across egg masses.
4. Treatment is tested within cohort first. Cross-cohort comparison is made on
   overlap and functional category, not a pooled model, since cohort is
   confounded with sequencing run.
5. A candidate novel protein family requires no database hit, support across
   multiple proteins and metagenomes, a length filter, removal of spurious-ORF
   matches, and exclusion of contig-edge fragments.
6. Language: cohort-associated or cohort-restricted protein space, never
   cohort-specific biology or year-driven expansion.
