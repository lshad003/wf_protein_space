# Wood frog gut protein space: analysis record

This repository documents the construction of a non-redundant protein catalog
from wood frog (*Lithobates sylvaticus*) faecal metagenomes collected across
three years, and the analyses used to characterize that protein space: its size
and saturation, the contribution of each cohort, the technical variables that
drive apparent richness, and its functional annotation.

## Sample set

89 faecal metagenomes from 64 animals. Faecal samples were collected weekly and
pooled by month; each metagenome is one monthly pool from one animal, with stems
named `<animal>.<timepoint>`.

| cohort | metagenomes | animals | structure |
|--------|-------------|---------|-----------|
| WF22 | 44 | 19 | multiple monthly pools per animal |
| WF23 | 9 | 9 | one monthly pool per animal, same month |
| WF24 | 36 | 36 | one monthly pool per animal |

WF22 came from three egg masses with two *Basidiobolus* strains against control
across three monthly timepoints; 6 of 50 metagenomes were excluded (one control
with contaminating *Basidiobolus* reads, one on quality review, four animals that
died before completing the experiment). WF23 sequenced animals are all egg mass
3, 5 control against 4 UHM520.7734. WF24 sequenced animals are 33 egg mass 4 and
3 egg mass 2, across 7 treatment codes with 5 animals each. Tank is confounded
with egg mass in WF24.

## Steps

### Step 1. Gene prediction and catalog input

Per-sample assembled contigs are translated with Prodigal in metagenomic mode,
and headers renamed so every protein carries its source metagenome. Predictions
existed for 31 of the 45 WF23 and WF24 metagenomes; the remaining 14 were
predicted with the same tool version and parameters.

**Main result.** 24,151,134 proteins from the 45 WF23 and WF24 metagenomes
(WF23 3,642,871; WF24 20,508,263), joining 16,399,461 from the 44 WF22
metagenomes for a catalog input of 40,550,595. Gene density per unit of
assembled sequence is constant across prediction batches (median 4,580 proteins
per Mb for pre-existing predictions against 4,602 for those predicted here,
range 4,047 to 4,969), so the two runs are interchangeable. UHM585.41009 is a
low outlier at 313,220 proteins from a 64 Mb assembly.

| File | Purpose |
|------|---------|
| `scripts/step0_check_existing_predictions.sh` | Which metagenomes already had predictions |
| `scripts/step0b_check_by_assembly.sh` | Per-stem inventory of existing protein files |
| `scripts/step0c_locate_predict_script.sh` | Locate prediction pipeline, confirm tool version |
| `scripts/step1b_preflight.sh` | Contig files confirmed for the 14 missing metagenomes |
| `scripts/step1c_prodigal_array.sh` | Prodigal, 14 metagenomes |
| `scripts/step1d_rename_45.sh` | Biosample-prefixed headers, all 45 |
| `scripts/step1e_protein_counts.sh` | Per-sample protein counts against assembly size |

Output: `results/step1_protein_counts_45.tsv`, `results/step0_by_assembly_check.tsv`

### Step 2. Catalog construction

All 40,550,595 proteins clustered at three amino-acid identity thresholds in one
run. A database completeness check runs before clustering, after a truncated
database was detected in an earlier attempt.

**Main result.** 6,182,117 representatives at 95% identity, 5,361,227 at 90%,
2,922,537 at 50%. Catalog is `catalog/db/LsPS_AA`.

| File | Purpose |
|------|---------|
| `scripts/step2_cluster.sh` | Input assembly, database build, clustering at three tiers |
| `scripts/step2a_preflight.sh` | Tool version, database format, disk headroom |

Output: `catalog/db/LsPS_AA_{95,90,50}_rep.fasta`

### Step 3. Search-based coverage of the earlier catalog

A random subsample of new-cohort proteins is searched against the earlier
WF22-only catalog to measure how much of the new protein space it already
represented. A positive control of WF22 proteins, which built that catalog and
must therefore be recovered, validates the search settings before the result is
interpreted.

**Main result.** Coverage mode determines the answer. Under target-side coverage the positive control recovered only 61.43% of proteins that are present in the catalog by construction, because short proteins cannot cover 80% of a longer representative; mapped proteins had median length 261 aa against 123 aa for unmapped. Under query-side coverage the control recovers 98.25%, and on that setting 46.15% of new-cohort proteins match the earlier catalog at gene level and 66.66% at family level.

| File | Purpose |
|---|---|
| `scripts/step2c_diagnostic_epyc.sh` | Subsample search, gene-level tier |
| `scripts/step2d_control_and_tiers.sh` | Positive control and family-level tier |
| `scripts/step2e_covmode2.sh` | Coverage-mode comparison with control |

### Step 4. Catalog composition and cluster support

Every cluster is resolved into the cohorts its members come from, at each tier,
before and after filters on cluster support. Singleton clusters are necessarily
cohort-exclusive, so cohort sharing cannot be interpreted until they are removed.

**Main result.** 55.27% of gene-level clusters are singletons. A support filter of at least three members from at least two metagenomes leaves 2,069,453 clusters, and cohort sharing rises sharply: clusters containing all three cohorts go from 6.71% to 20.06%, and at family level to 26.00%. Cluster size is heavily skewed, with the top 1% of family-level clusters holding 45.85% of all proteins.

| File | Purpose |
|---|---|
| `scripts/step3_cluster_composition.sh` | Cluster membership tables and cohort composition |
| `scripts/step4_nonsingleton.sh` | Composition under four support filters, size distribution |

### Step 5. Cohort contribution at matched sampling effort

Cohorts differ in metagenome count, so cohort-exclusive cluster counts partly
measure sequencing effort rather than cohort identity. Cohorts are compared at
equal numbers of metagenomes across repeated random draws.

**Main result.** At 9 metagenomes per cohort across 10 draws, WF22-exclusive clusters fall from 15.12% to 9.26% once the sample-count advantage is removed, while WF24-exclusive remains 26.13% and clusters shared by all three cohorts are 24.22%. Per-cohort richness at 9 metagenomes is 765,122 clusters for WF22, 1,166,662 for WF23 and 2,331,079 for WF24, tracking per-sample protein yield rather than cohort identity.

| File | Purpose |
|---|---|
| `scripts/step5_equaln_rarefaction.sh` | Equal-n cohort comparison and per-cohort accumulation |

### Step 6. Batch structure

Sequencing run, centre, platform and read depth are tabulated for all 89
metagenomes, so that cohort comparisons can be checked against the batch
structure they are confounded with.

**Main result.** Sequencing runs are not fully confounded with cohort. WF23 and WF24 share one run, and WF24 spans two runs, which permits direct estimation of batch effects. Within the shared run WF23 and WF24 yield 422,857 and 415,208 proteins per metagenome, so cohort alone does not affect protein yield. Across its two runs WF24 yields 415,208 against 756,692 proteins per metagenome, a 1.8-fold difference on 1.24-fold more reads. The difference is assembly size rather than gene prediction, since gene density per Mb is constant across runs. WF22 occupies separate runs at two centres plus one AVITI run, so WF22 against new-cohort contrasts cannot be separated from batch.

| File | Purpose |
|---|---|
| `scripts/step6a_inspect_qc.sh` | Read statistics format and run assignment |
| `scripts/step6b_batch_table.sh` | Batch table and the two key contrasts |

Output: `results/batch_table_89.tsv`

### Step 7. Depth-normalized rarefaction

Because assembled sequence per metagenome is the dominant technical variable,
accumulation is measured against proteins sampled rather than metagenomes
sampled, and WF24 is split by sequencing run so that cohort and run are not
conflated.

**Main result.** At 3,600,000 proteins WF22 gives 1,169,677 clusters and WF23 1,159,262, a 0.9% difference. WF24 gives 1,285,504 in one run and 1,687,560 in the other, so the same cohort differs by 31% across runs while different cohorts within one run differ by 11%. Sequencing run, not cohort, is the dominant driver of catalog richness. The pooled curve reaches 6,182,117 clusters at 40,550,595 proteins and is still adding roughly 435,000 clusters per 5 million proteins, so gene-level protein space is not saturated.

| File | Purpose |
|---|---|
| `scripts/step7_rarefaction.sh` | Accumulation against proteins sampled, per cohort and pooled |
| `scripts/step7b_wf24_byrun.sh` | Richness at matched depth, WF24 split by run |

Output: `results/rarefaction_95.tsv`

### Step 8. Functional annotation

Gene-level representatives are annotated against eggNOG to establish the
annotated and unannotated fractions of the catalog. The reference database is
staged to node-local storage, since concurrent access to the shared copy stalls
the annotation phase entirely. Annotation from a stored hits table fails in both
available versions of the mapper, so every chunk is run end to end.

| File | Purpose |
|---|---|
| `scripts/step8a_annot_preflight.sh` | Available tools, databases and reference invocations |
| `scripts/step8e_test.sh` | Full run against hits-table reuse, two mapper versions |
| `scripts/step8f_eggnog_v2.sh` | eggNOG annotation, database staged node-local |
| `scripts/step8c_pfam_array.sh` | Pfam annotation |

## Planned steps

| Step | Description |
|---|---|
| 9 | Known and unknown fractions, unknowns split by whether they match sequenced genomes |
| 10 | Quality control of unknowns: length, spurious-ORF screen, complete against partial calls |
| 11 | Novel family calling on supported clusters with no database match |
| 12 | Structure prediction and structure search for novel family representatives |
| 13 | Taxonomic composition against read profiles |
| 14 | Abundance across all 89 metagenomes; prevalence, core and accessory, treatment contrasts |

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

## INSTRUCTIONS FOR CLAUDE CODE
Execute one step per session, the one named. After its checks pass, append a
dated entry to logs/step_log.md and stop. Never modify anything under
/bigdata/stajichlab/shared/. Write outputs under this project directory. Heavy
work goes in an sbatch script on -p epyc, never the login node. Do not commit
data files. Before stating any number, show the command that produced it.
