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

Proteins were predicted from assembled contigs with Prodigal in metagenomic mode.

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

### Step 11. Taxonomic composition of the annotated fraction

Broad clade is read from the eggNOG orthologous group assignments of the
annotated representatives, and examined against cluster size, since host and
dietary sequence would be expected to concentrate in poorly supported clusters.

**Main result.** The annotated fraction is overwhelmingly bacterial: 92.64% Bacteria, 6.93% Eukaryota, 0.27% Archaea and 0.16% Viruses of 2,980,052 annotated representatives. The most frequent assignment levels are Alphaproteobacteria 16.30%, Bacteroidetes 13.95%, Actinobacteria 13.34%, Gammaproteobacteria 11.44% and Betaproteobacteria 10.31%, with Metazoa at 3.05% and Fungi at 2.57%. Eukaryotic assignment declines with cluster support, from 7.63% among singletons to 3.61% among clusters of twenty-one or more, consistent with host and dietary sequence being present but sparse and poorly replicated. This applies to annotated proteins only; the unannotated fraction carries no taxonomic assignment and its composition is not established by this analysis.

| File | Purpose |
|---|---|
| `scripts/step13_taxonomy.sh` | Broad clade and assignment level, against cluster size; 89 metagenomes, superseded by `scripts/step13_taxonomy_n88.sh` |
| `scripts/step13_taxonomy_n88.sh` | Broad clade and assignment level, 88 metagenomes |

Output: `results/taxonomy_summary.tsv`; at 88: `results/taxonomy_summary_n88.tsv`

### Step 12. Annotation coverage and the four-way classification

Two annotation sources are combined before any dark fraction is quoted, since
eggNOG alone overstates it. Pfam is assigned with hmmsearch under gathering
thresholds, which fixes the criterion without a chosen E-value. The union then
defines what counts as known, and everything unknown is carried forward.

**Coverage.** Of 6,171,602 representatives, eggNOG annotates 2,980,052
(48.29%), Pfam 2,709,544 (43.90%), and their union 3,222,677 (52.22%). Within
the support-filtered set of 2,067,011, the union leaves 585,347 unannotated
(28.32%). At the 50% family tier the supported dark fraction falls from 55.40%
under eggNOG alone to 51.47% under the union, so about four points of the
earlier figure were an artefact of using one database.

**Classes.** Representatives are labelled K if they carry a Pfam domain,
KWP if eggNOG assigns them but Pfam does not, and U if neither does. Supported
counts are K 1,343,020, KWP 138,644, U 585,347. Of the K representatives,
76,515 carry only DUF or UPF domains, which a stricter definition would move
out of K, leaving 1,266,505.

| File | Purpose |
|---|---|
| scripts/step12c_check_pfam_done.sh | Pfam completion verified before any Pfam number is used |
| scripts/step12d_union_replevel.sh | Representative-level union of eggNOG and Pfam; 89 metagenomes, superseded by scripts/step12d_union_replevel_n88.sh |
| scripts/step12d_union_replevel_n88.sh | Representative-level union, 88 metagenomes |
| scripts/step12e_family_union.py | Union carried to the 50% family tier; 89 metagenomes, superseded by scripts/step12e_family_union_n88.py |
| scripts/step12e_family_union_n88.py | Union carried to the 50% family tier, 88 metagenomes |
| scripts/step12f_rep_classes.sh | Per-representative K, KWP and U labels with the DUF split; 89 metagenomes, superseded by scripts/step12f_rep_classes_n88.sh |
| scripts/step12f_rep_classes_n88.sh | K, KWP and U labels, 88 metagenomes; classes checked equal to the 89 labels |

Output: results/union_rep_level.tsv, results/family_union_dark.tsv,
results/rep_classes.tsv, results/known_union_ids.txt; at 88:
results/union_rep_level_n88.tsv, results/family_union_dark_n88.tsv,
results/rep_classes_n88.tsv, results/known_union_ids_n88.txt

### Step 13. Genomic and environmental unknowns

The unannotated representatives are searched against NCBI ClusteredNR to
separate genes that exist in sequenced genomes but have no assigned function
from genes that appear nowhere. The database is hard-coded after inspecting it
with dbinfo rather than selected by a glob, because the database defines the
boundary being reported. Hits are recorded loosely and the criterion is applied
afterwards, so alternatives can be compared without repeating the search.

**The search.** 586,215 supported unknown representatives against
nr_cluster_seq 20260128 (470,748,714 sequences) with DIAMOND 2.1.24
--very-sensitive, recording up to five targets per query at E <= 1e-3.
Thirteen hours on 32 cores with the database staged to node-local scratch.
At 88 metagenomes 585,347 of these remain supported and are classified from
the same search.

**The criterion matters.** Loose to strict, the environmental unknown count
runs 374,133 (63.92%) at E <= 1e-5 with no coverage requirement, 430,301
(73.51%) at E <= 1e-10 with query and subject coverage at least 50%, and
448,567 (76.63%) at E <= 1e-20 under the same coverage. The middle criterion
is adopted, since coverage on both sides prevents a short shared motif from
counting as a match.

**Main result.** Of the supported catalogue, 430,299 representatives (20.82%)
have no hit in a database of 470 million proteins. The genomic unknowns,
155,048 representatives (7.50%), mostly match sequences that are themselves
labelled hypothetical, so they are recognised without being characterised.
Mobile element titles account for 3,954 of the 211,712 loose hits (1.9%).

| File | Purpose |
|---|---|
| scripts/step13e_db_readability.sh | Candidate databases inventoried and checked for readability |
| scripts/step13d_extract_u.py | Unannotated representatives extracted to fasta |
| scripts/step13f_nrclust_search.sh | Search against ClusteredNR submitted |
| scripts/step13h_gu_eu_split.sh | Split reported under several stated criteria; 89 metagenomes, superseded by scripts/step13h_gu_eu_split_n88.sh |
| scripts/step13h_gu_eu_split_n88.sh | Split under the same criteria, supported unknowns at 88 metagenomes |
| scripts/step13h_gu_eu_split_89log.sh | Unchanged step13h run once more to record its 89-sample output in a log |
| scripts/step13i_fourway.py | Four-way classification assembled with per-class properties; 89 metagenomes, superseded by scripts/step13i_fourway_n88.py |
| scripts/step13i_fourway_n88.py | Four-way classification, 88 metagenomes |

Output: results/step13_nrclust20260128_hits.tsv, results/fourway_classes.tsv; at 88:
results/fourway_classes_n88.tsv

### Step 14. Structure of the unknown fraction

Two controls are applied to the unannotated set before it is described as
biology. AntiFam screens for spurious open reading frames using the same
profiles and version used by the global survey this work is compared against.
Clustering then asks whether unknown proteins resemble each other, which is the
step that would produce protein families if they existed here.

**Spurious sequences are rare.** AntiFam v6.0 flags 1,464 of 585,347
representatives, 0.2501%. The published global survey flagged 43 of 19,986,348,
but screened a set already filtered to families of 100 or more members, so the
two rates are not directly comparable and both criteria are stated.

**Unknown proteins do not form families.** Clustered at 50% identity with 80%
coverage, the 585,347 unknown representatives give 536,030 families. Only 8,464
have three or more members and 10 have 100 or more. At 88 metagenomes these are
the families of the 89-sample clustering restricted to supported members, not a
new clustering. The same measure applied to
the 1,344,403 known representatives gives 398,126 families, 102,134 with three
or more members and 379 with 100 or more. Known genes therefore cluster about
sixteen times more often on identical data with an identical method, so the
result is a property of the unknown fraction rather than of sample size.

**Consequence for framing.** A family catalogue comparable to the 106,198
novel metagenome protein families reported from 26,931 metagenomes cannot be
built from 88. What this catalogue documents is non-redundant unknown protein
space in one host, not a set of novel families.

| File | Purpose |
|---|---|
| scripts/step14b_antifam.sh | AntiFam screen of the unannotated representatives; 89 metagenomes, superseded by scripts/step14b_antifam_n88.sh |
| scripts/step14b_antifam_n88.sh | AntiFam screen rerun on the 585,347 supported unknowns at 88 metagenomes |
| scripts/step14c_novel_families.sh | Unknown representatives clustered at 50% identity; 89 metagenomes, superseded by scripts/step14c_novel_families_n88.sh |
| scripts/step14c_novel_families_n88.sh | Family sizes at 88 metagenomes from the 50% clustering, no reclustering |
| scripts/step14d_unk_40.sh | Clustering repeated at 40% identity as a threshold control |

Output: results/antifam_hits.tblout, results/unk_clusters_50.tsv,
results/unk_family_sizes.tsv; at 88: results/antifam_hits_n88.tblout,
results/unk_clusters_50_n88.tsv, results/unk_family_sizes_n88.tsv

### Step 16. Abundance across all 88 metagenomes

Reads are mapped to the support-filtered catalogue rather than the full one,
since singleton representatives are contig-edge fragments and would compete for
reads belonging to the intact gene. Two counting criteria are written for every
gene in the same pass, so the choice between them is made at analysis time and
never requires a remap. No alignment file is written to disk.

**Reference.** The 2,069,453 supported representatives as nucleotide coding
sequences, 1,750,474,248 bp, mean 845.9 bp, indexed with bwa-mem2 2.3. Reads
are resolved through the sequencing manifest, with all 89 samples matched to
distinct read pairs.

**Mapping is even across cohorts.** Mean mapped fraction 0.9587 in WF22
(n = 44), 0.9634 in WF23 (n = 9) and 0.9667 in WF24 (n = 35), a spread under
one percentage point. The catalogue therefore represents all three collection
years equally well, and the unmapped remainder is a measure of how much of each
metagenome falls outside the analytical catalogue.

**Verification.** Column sums of both count matrices reproduce the per-sample
mapped-read totals exactly for all 89 samples. Genes detected at least once:
2,066,746 under primary counting and 2,040,274 under the MAPQ >= 10 subset.

| File | Purpose |
|---|---|
| scripts/step16b_supported_set.py | Support-filtered representative set materialized and cross-checked |
| scripts/step16e_wf22_cds_rename.sh | WF22 coding sequences renamed and verified against the catalogue identifiers |
| scripts/step16h_supported_cds.py | Mapping reference extracted for the supported representatives |
| scripts/step16j_index_reference.sh | Reference length table built and indexed |
| scripts/step16k_manifest_join.sh | Read files resolved for all 89 samples through the manifest |
| scripts/step16q_map_array.sh | Mapping array, both counting criteria, no alignment file retained |
| scripts/step16t_qc_table.sh | Per-sample mapping quality table with cohort labels; 89 metagenomes, superseded by scripts/step16t_qc_table_n88.sh |
| scripts/step16t_qc_table_n88.sh | Mapping quality table without UHM586.41010, by cohort |
| scripts/step16u_matrix.py | Count matrices assembled and checked against the summaries; 89 metagenomes, superseded by scripts/step16u_matrix_n88.py |
| scripts/step16u_matrix_n88.py | Both matrices at 88 metagenomes from the 89 matrices; column sums checked |

Output: results/count_matrix_primary.tsv.gz, results/count_matrix_mapq10.tsv.gz,
results/mapping_qc_89.tsv, results/supported_reps_95.txt

At 88 metagenomes the analyses use results/count_matrix_primary_n88.tsv.gz:
the UHM586.41010 column removed and rows restricted to the 2,067,011 clusters
still supported (scripts/step25a_filter_n88.py). Reads were not remapped, so
the reference and the column-sum verification above describe the 89-sample
mapping; the mapped fractions and detected-gene counts are given for the 88
metagenomes (results/mapping_qc_88.tsv, results/count_matrix_mapq10_n88.tsv.gz).

### Step 17. Treatment, maternal origin and time

The WF22 design is longitudinal: 44 samples come from 15 animals sampled across
three months. Any model treating those samples as independent inflates its
p-values, so animal is used as a blocking factor. Each cohort is tested
separately, since cohort is confounded with sequencing run. Genes are filtered
to those present in most samples of the cohort before testing.

**Method.** limma-voom with TMM normalization. Repeated sampling in WF22 is
handled with duplicateCorrelation blocking on animal. LinDA was attempted first
and abandoned: it allocates a matrix quadratic in the number of features and
requested 1,010 Gb for 368,236 genes.

**Ignoring the repeated sampling changes the answer.** The consensus
within-animal correlation is 0.2864. Without blocking, the WF22 treatment
contrasts return 5,005 and 7,263 genes. With blocking, on identical data, they
return 3 and 155.

**No treatment effect in any cohort.** Within egg mass 3, where the design is
balanced, STP1710.7 against control returns zero genes from 21 samples and
seven animals. WF23 returns zero for UHM520.7734 against control. WF24 returns
zero for five of six strains and three genes for the sixth (code 3, UHM516.7697; 34 animals). Codes 5 and 6 in WF24
are the two strains used in WF22, so the null is replicated in a second year
with different animals. Across nine contrasts, three cohorts and seven
Basidiobolus strains, gene abundance does not respond detectably.

**Maternal origin and time do.** In the same blocked WF22 model, egg mass
returns 57,120 and 51,808 genes and month returns up to 141,580. The negative
treatment result is therefore not a failure of power in the design.

**An aggregate test agrees.** The unannotated share of gene abundance does not
differ by treatment in any cohort (Kruskal-Wallis p = 0.4688 and 0.4222 in WF22,
Mann-Whitney p = 0.1111 and 0.1905 in WF23, Kruskal-Wallis p = 0.4481 and 0.5527
in WF24 with 34 animals).

| File | Purpose |
|---|---|
| scripts/step17a_sample_table.sh | Treatment and egg mass assembled for all 89 samples |
| scripts/step17b_class_shift.py | Unannotated abundance share tested by treatment within cohort; 89 metagenomes, superseded by scripts/step17b_class_shift_n88.py |
| scripts/step17b_class_shift_n88.py | Unannotated abundance share by treatment, 88 metagenomes, WF24 34 animals |
| scripts/step17b_class_shift_89log.sh | Unchanged step17b run once more to record its 89-sample output in a log |
| scripts/step17c_wf22_meta.sh | WF22 design table with animal, month and egg mass |
| scripts/step17g_limma_wf22.R | WF22 tested blocked and unblocked, and within egg mass 3; 89 metagenomes, superseded by scripts/step17g_limma_wf22_n88.R |
| scripts/step17g_limma_wf22_n88.R | WF22 tested blocked and unblocked, and within egg mass 3, 88-metagenome matrix |
| scripts/step17h_limma_wf2324.R | WF23 treatment contrasts; its WF24 run (36 animals) is superseded by step17i; 89 metagenomes, superseded by scripts/step17h_limma_wf2324_n88.R |
| scripts/step17h_limma_wf2324_n88.R | WF23 treatment contrasts, 88-metagenome matrix |
| scripts/step17i_limma_wf24_n34.R | WF24 treatment contrasts, 34 animals, exclusions in metadata/wf24_excluded.tsv; 89 metagenomes, superseded by scripts/step17i_limma_wf24_n34_n88.R |
| scripts/step17i_limma_wf24_n34_n88.R | WF24 treatment contrasts, 34 animals, 88-metagenome matrix |

Output: results/limma_A_blocked_*.csv, results/limma_B_naive_*.csv,
results/limma_C_EM3_*.csv, results/limma_WF23_*.csv, results/limma_WF24_n34_*.csv (results/limma_WF24_treatment*.csv
are the superseded 36-animal run),
metadata/wf22_design.tsv, metadata/wf24_treatment_key.tsv, metadata/wf24_excluded.tsv;
at 88: results/limma_*_n88.csv, metadata/samples_88_treatment.tsv

### Step 19. Properties of the unknown fraction

Whether the unannotated genes are biology or artefact is settled with measured
properties rather than argument. Length, open reading frame completeness,
prevalence, persistence across collection years and abundance share are compared
across the four classes on the same catalogue.

**They are short, and not because they are broken.** Mean coding length is
1,045 bp for K, 795 for KWP, 498 for GU and 367 for EU. Restricting to complete
open reading frames barely moves the unannotated figure, from 402 to 405 bp,
while the known figure rises from 1,045 to 1,093. Complete reading frames make
up 79.5% of K, 69.8% of KWP, 76.0% of GU and 63.7% of EU. Unannotated genes are
therefore small complete proteins of about 134 amino acids, a size class that
reference databases represent poorly.

**They are widespread.** 82.77% of unannotated genes occur in more than 20 of
88 samples and 3,496 occur in all 88. Only 0.65% occur in five or fewer.

**They persist across years.** 496,266 unannotated genes occur in all three
collection years and carry 91.15% of unannotated gene reads. The 18,896
single-year genes are shorter, mean 246 bp, and carry 1.23%. A further 260
genes attract no reads at all and are assembly artefacts. Requiring presence in
all three years therefore removes almost no signal, which makes that subset a
defensible core set for claims about persistence.

**Abundance depends on how it is measured.** Unannotated genes are 28.32% of
supported genes but 7.04% of raw mapped reads. After length normalization they
are 19.32% of gene copies. The gap is entirely explained by their short length,
so the length-normalized figure is the one that describes the community and the
raw figure is reported alongside it. The share is stable across cohorts, 7.11%
in WF22, 6.28% in WF23 and 7.16% in WF24, with a per-sample range of 4.38 to
16.12%.

| File | Purpose |
|---|---|
| scripts/step16v_prevalence.py | Prevalence per gene by cohort and class; 89 metagenomes, superseded by scripts/step16v_prevalence_n88.py |
| scripts/step16v_prevalence_n88.py | Prevalence per gene by cohort and class, 88 metagenomes |
| scripts/step16w_dark_abundance.py | Abundance share of each class, raw and length-normalized; 89 metagenomes, superseded by scripts/step16w_dark_abundance_n88.py |
| scripts/step16w_dark_abundance_n88.py | Abundance share of each class, 88 metagenomes |
| scripts/step16y_complete_by_class.py | Reading frame completeness by class from the gene caller flags; 89 metagenomes, superseded by scripts/step16y_complete_by_class_n88.py |
| scripts/step16y_complete_by_class_n88.py | Completeness and mean length by class, 88 metagenomes |
| scripts/step18a_dump_completeness.py | Completeness flags written per representative |
| scripts/step19a_core_unknown.py | Unannotated genes characterized by number of years detected; 89 metagenomes, superseded by scripts/step19a_core_unknown_n88.py |
| scripts/step19a_core_unknown_n88.py | Unannotated genes by number of years detected, 88 metagenomes |
| scripts/step18b_fig12.py | Prevalence and abundance figures |

Output: results/prevalence_primary.tsv.gz, results/dark_abundance_by_sample.tsv,
results/completeness_by_class.tsv, results/core_unknown_by_years.tsv.gz,
figures fig_prevalence_by_class.png and fig_dark_abundance.png; at 88:
results/prevalence_primary_n88.tsv.gz, results/dark_abundance_by_sample_n88.tsv,
results/core_unknown_by_years_n88.tsv.gz, results/completeness_by_class_n88.tsv

### Step 20. Taxonomic context of the unknown fraction

Short unannotated proteins carry little taxonomic signal on their own, so
taxonomy is inherited from the contig each gene sits on, using assignments made
upstream against UniRef50. The question is not what species these genes come
from, which the data cannot answer, but whether unannotated genes sit in more
taxonomically obscure genomic neighbourhoods than annotated ones.

**Coverage.** Contig assignments exist for 48 of 88 samples, giving taxonomy
for 888,183 of 2,067,011 genes, 42.97%. Every figure below is conditional on
that subset and the comparison between classes carries the result, not the
absolute values.

**Unannotated genes sit on more obscure contigs.** Contigs are unclassified for
74.6% of K genes, 74.6% of KWP, 74.4% of GU and 91.6% of EU. Among the genes
whose contigs are classified, the environmental unknowns also invert the usual
pattern: 4.5% eukaryotic against 3.9% bacterial, where the other three classes
run about 23% bacterial and 2% eukaryotic. Genes with no database match are
therefore embedded in DNA that is itself taxonomically orphan, which is what
lineages without sequenced relatives would produce.

| File | Purpose |
|---|---|
| scripts/step20a_tax_recon.sh | Existing scaffold classifications located and checked |
| scripts/step20b_contig_tax.sh | Assignment file format and coverage verified |
| scripts/step20c_gene_tax.py | Contig taxonomy inherited by gene and summarized by class; 89 metagenomes, superseded by scripts/step20c_gene_tax_n88.py |
| scripts/step20c_gene_tax_n88.py | Contig taxonomy by gene, 88 metagenomes; reads .tsv or .tsv.gz contig files |

Output: results/gene_taxonomy.tsv; at 88: results/gene_taxonomy_n88.tsv

### Step 21. Gene neighbourhood of unknown families

Genes numbered consecutively on a contig are physical neighbours, and bacterial
genes in the same pathway tend to sit together. Where an unknown family
repeatedly appears beside the same annotated domain, that neighbour is a
functional hint. This is applied only to the 8,464 unknown families with three
or more members, since a family seen on one contig cannot show a conserved
neighbourhood. A member counts once for a domain regardless of whether it
appears on one side or both.

**Nearly half have a conserved neighbour.** Of the 8,464 testable families,
3,614 (42.7%) have an annotated domain beside at least half their members. The
strongest cases, families of eight to ten members with a neighbour present in
every one, sit beside oxidoreduction and respiration domains (GSDH, Rieske,
COX1, adh_short, Aldo_ket_red), transporters (MFS_1, ATP_bind_1, mechanosensitive
channels), regulators (HTH_1, TetR-like, PhoU), outer membrane and surface
proteins (Fimbrial, AsmA), and peptidases. The unannotated genes are therefore
embedded in ordinary bacterial metabolic and regulatory machinery rather than in
mobile elements.

**Two limits.** This applies only to families with three or more members, whose
members are 6.29% of supported unannotated genes, so it describes the testable minority. And
adjacency is read from gene numbering on a contig, so it inherits any assembly
error in that contig. A conserved neighbour is a hypothesis about function, not
an assignment.

| File | Purpose |
|---|---|
| scripts/step21a_neighbour_check.sh | Gene identifiers confirmed to encode contig and position; 89 metagenomes, superseded by scripts/step21a_neighbour_check_n88.sh |
| scripts/step21a_neighbour_check_n88.sh | Identifier check on the 88-metagenome tables |
| scripts/step21b_neighbours.py | Annotated domains counted beside members of each unknown family; 89 metagenomes, superseded by scripts/step21b_neighbours_n88.py |
| scripts/step21b_neighbours_n88.py | Neighbour domains of unknown families, 88 metagenomes |

Output: results/unk_family_neighbours.tsv; at 88: results/unk_family_neighbours_n88.tsv

### Step 22. Community-level and targeted tests of treatment

Per-gene tests answer whether individual genes respond. They do not answer
whether the community as a whole shifts, and a signal confined to a small
functional class could be diluted across hundreds of thousands of genes. Three
further designs address both, and the permutation scheme is treated as part of
the design rather than a default.

**Permutation has to match the design.** Treatment and egg mass are constant
within an animal, so a scheme that permutes freely inside animal blocks never
moves those labels and returns a meaningless result. Between-animal terms are
therefore tested by permuting whole animals, which vegan allows only on a
balanced design, so the one animal with two samples is dropped, leaving 42
samples from 14 animals. Month varies within an animal and is permuted inside
animals. An earlier run using block permutation for all three terms is
superseded.

**Composition does not shift with treatment.** On Bray-Curtis distances over
Hellinger-transformed abundances, treatment explains 5.1% of variation with
whole-animal permutation (p = 0.928) and 6.9% on animal centroids, which
collapses the repeated measures to 15 independent units (p = 0.657). In WF24,
where 34 animals are each sampled once, treatment R2 is 0.173 (p = 0.679), at
the chance level expected for 6 of 33 degrees of freedom.
Constrained ordination on treatment, conditioning out egg mass and month,
recovers no significant axis.

**Nor does gene presence.** Scoring a gene as present at five or more reads
gives 1,637,544 genes that vary across samples, a mean of 523,226 per sample.
On Jaccard distances, treatment explains 6.2% (p = 0.665) and 10.6% on animal
centroids (p = 0.465), and in WF24 R2 0.170 (p = 0.770, 34 animals, chance level). The number of genes
detected per sample does not differ by treatment (p = 0.3676).

**Nor the carbohydrate subset, where an effect was most plausible.** No CAZy
annotation exists for this catalogue, so the subset is defined by Pfam family
name, a looser proxy: 59,297 representatives, of which 7,891 are present in at
least 35 of 44 samples. Blocked per-gene testing returns zero genes for both
treatments, against 1,541 and 2,046 for egg mass and 1,761 and 2,586 for month
in the same subset. Composition of the subset gives treatment 5.3%
(p = 0.454) against egg mass 67.3% (p = 0.001), and total carbohydrate
abundance does not differ (p = 0.3181).

**Maternal origin does shift composition.** Egg mass explains 32.7% of
variation under whole-animal permutation and 51.3% on animal centroids, both
p = 0.001, and its group dispersions do not differ (p = 0.214), so this is a
difference in composition rather than in variability. Month explains 12.7%
(p = 0.001).

**One caveat on the negative.** Treatment groups do differ in dispersion
(p = 0.003) while their centroids do not, so the groups vary in how much they
vary. This is reported rather than interpreted.

| File | Purpose |
|---|---|
| scripts/step22b_permanova_fix.R | Composition tested under three permutation designs; 89 metagenomes, superseded by scripts/step22b_permanova_fix_n88.R |
| scripts/step22b_permanova_fix_n88.R | Composition under three permutation designs, 88-metagenome matrix |
| scripts/step22c_jaccard.R | Presence/absence composition and gene richness; 89 metagenomes, superseded by scripts/step22c_jaccard_n88.R |
| scripts/step22c_jaccard_n88.R | Presence/absence composition and gene richness, 88-metagenome matrix |
| scripts/step22e_permanova_wf24_n34.R | WF24 abundance and presence/absence composition, 34 animals; 89 metagenomes, superseded by scripts/step22e_permanova_wf24_n34_n88.R |
| scripts/step22e_permanova_wf24_n34_n88.R | WF24 composition, 34 animals, 88-metagenome matrix |
| scripts/step22d_cazy_subset.R | Carbohydrate subset tested per gene and by composition; 89 metagenomes, superseded by scripts/step22d_cazy_subset_n88.R |
| scripts/step22d_cazy_subset_n88.R | Carbohydrate subset per gene and by composition, 88-metagenome matrix |

Output: logs step22b, step22c, step22d and step22e; results/limma_cazy_*.csv;
at 88: logs step22b_n88 to step22e_n88, results/limma_cazy_*_n88.csv

## Planned steps

| Step | Description |
|---|---|
| 12 | Pfam domain assignment and the union of eggNOG and Pfam coverage |
| 13 | Four-way split: known with domain, known without domain, genomic unknown, environmental unknown |
| 14 | Novel family calling on supported clusters with no database match |
| 15 | Structure prediction and structure search for novel family representatives |
| 16 | Abundance across all 88 metagenomes; prevalence, core and accessory |
| 17 | Treatment contrasts within cohort |

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

No count is stated without the output file that produced it. MMseqs2 is pinned
to 13-45111 and runs only on the epyc partition, since other nodes fail with an
illegal instruction. Coverage mode is stated explicitly for every search, since
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
