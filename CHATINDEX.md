# chat_index.md - WF Protein Space
# APPEND ONLY after 2026-08-21. New entries at the bottom via cat >>.
# Entry format:
# ## [CATEGORY] title (date)
#   - key results / decisions
#   - files produced (full paths)

## [ANALYSIS] WF23/WF24 integration - scope and catalog strategy (Aug 2026)
- Settled framework: map-to-existing diagnostic first; rebuild costs all verified IDs;
  clusterupdate is the middle path (pin MMseqs2 v13).
- Located year databases in /bigdata/stajichlab/lshad003/ncbi-deposit/wf_databases/

## [CODE] sample manifest - wf_year counts and stem lists (Aug 21, 2026)
- Master list: /bigdata/stajichlab/lshad003/ncbi-deposit/results/deposit_classification_v2.tsv
- Year column is wf_year; classification is only WF-vs-wild (bug hit twice, documented)
- Verified: WF22 50/19 longitudinal; WF23 9/9; WF24 36/36; wild 79 NA
- Wrote wf22/23/24_stems.txt under ncbi-deposit/results/, copied to
  /bigdata/stajichlab/lshad003/wf_protein_space/metadata/

## [ANALYSIS] WF23/24 metadata join - design coverage of sequenced samples (Aug 21, 2026)
- CORRECTION: "single timepoint" describes SEQUENCED data only; colonies sampled
  longitudinally (WF23 median 38 dates/animal; WF24 weekly pools), banked in freezer.
- Joined all 45 sequenced stems to year databases, 45/45 matched:
  /bigdata/stajichlab/lshad003/wf_protein_space/metadata/wf23_wf24_sequenced_metadata.tsv
- CORRECTION: sequenced WF23 all egg_mass 3; WF24 33/36 egg_mass 4.
  Sequenced new years support catalog growth + treatment only; no maternal axis.
- BLOCKER: library number (40991+) to lab tube number (7870s/8270s) map missing;
  collection dates unknown for all 45; ask Walker group for submission sheet.
- Drive MIMS table confirmed stem-to-year for all 45, collection_date blank.

## NEXT CHATS TO OPEN
- [CODE] Prodigal array - predict proteins for 45 WF23/WF24 stems
- [ANALYSIS] LsFMGC95 map-to-existing diagnostic - WF23/24 mapping rate

## [ANALYSIS] Catalog build through annotation (Aug 21-23, 2026)
- Rebuilt the catalog from scratch as one three-cohort resource: 89 metagenomes
  (44 WF22 + 9 WF23 + 36 WF24), 40,550,595 input proteins.
  Tiers: 6,182,117 at 95%, 5,361,227 at 90%, 2,922,537 at 50%. 100% tier skipped.
- Predicted proteins for the 14 WF24 stems that lacked them (Prodigal v2.6.3),
  renamed all 45 new stems to biosample-prefixed headers.
- Search diagnostics: --cov-mode 1 fails on short proteins (WF22 positive control
  recovers only 61.43%); --cov-mode 2 recovers 98.25% and is the correct setting.
- Composition: 55.27% of gene-level clusters are singletons; support filtering
  (>=3 members, >=2 metagenomes) leaves 2,069,453 clusters and raises three-cohort
  sharing from 6.71% to 20.06%.
- Batch: sequencing run outweighs cohort. Same cohort across two runs differs 31%
  at matched depth; different cohorts within one run differ 11%. Cause is assembly
  size, not gene prediction (gene density constant at ~4,600 proteins/Mb).
- ORF completeness: 46.91% of all proteins complete; representatives 35.28%;
  rises from 11.48% in singletons to 86.14% at 21+ members, so the singleton
  fraction is contig-edge fragmentation.
- eggNOG: 48.27% of representatives annotated (v1 WF22-only was 47.2%).
  Family tier dark fraction 76.01% raw, 55.35% after support filtering.
  Annotation does NOT track ORF completeness: proteins spanning short contigs
  annotate better than complete ORFs within every cluster-size bin below 21.
- Taxonomy of the annotated fraction: 92.64% Bacteria, 6.93% Eukaryota,
  0.27% Archaea, 0.16% Viruses. Eukaryotic share falls from 7.64% in singletons
  to 3.59% at 21+ members. Fungi 2.57% (76,670 reps) worth checking against
  Basidiobolus treatment later.
- Repo: github.com/lshad003/wf_protein_space, README restructured as a
  step-by-step analysis record (Steps 1-11), gitignore now a whitelist.
- PENDING at close: Pfam via hmmsearch, jobs 27712560 and 27712562,
  output results/pfam_hs, ~15h remaining.
- NEXT: eggNOG+Pfam union, AGNOSTOS four-way split, then abundance on the
  support-filtered catalog and within-cohort treatment contrasts.

# ==============================================================
# From here on this file is CHATINDEX.md: append-only current state.
# Newest block wins. A claim named in a SUPERSEDES line is dead.
# Everything above is the earlier chat_index and is preserved as history.
# ==============================================================

## 2026-08-23 catalog build through annotation

VERIFIED (file path behind each number)
- Catalog: 89 metagenomes (44 WF22 + 9 WF23 + 36 WF24), 40,550,595 input proteins.
  Tiers 6,182,117 at 95%, 5,361,227 at 90%, 2,922,537 at 50%.
  catalog/db/LsPS_AA_{95,90,50}_rep.fasta
- Gene density constant across prediction batches, about 4,600 proteins per Mb.
  results/step1_protein_counts_45.tsv
- Searches require --cov-mode 2. WF22 positive control recovers 98.25% against
  61.43% under --cov-mode 1. logs/step2e.*.log
- 55.27% of gene-level clusters are singletons. Support filter (>=3 members from
  >=2 metagenomes) leaves 2,069,453 clusters and raises three-cohort sharing from
  6.71% to 20.06%. logs/step4.*.log
- Sequencing run outweighs cohort: same cohort across two runs differs 31% at
  matched depth, different cohorts within one run differ 11%.
  results/batch_table_89.tsv, results/rarefaction_95.tsv
- ORF completeness 46.91% overall, 35.28% of representatives, rising from 11.48%
  in singletons to 86.14% at 21+ members. results/orf_completeness.tsv
- eggNOG 48.27% of representatives annotated. Family-tier dark fraction 76.01%
  raw, 55.35% after support filtering across 1,022,162 families.
  results/annotation_summary.tsv, results/family_dark_fraction.tsv
- Annotated fraction is 92.64% Bacteria, 6.93% Eukaryota, 0.27% Archaea,
  0.16% Viruses. results/taxonomy_summary.tsv

SUPERSEDES
- Dead claim: "WF23 and WF24 are single timepoint per animal." Animals were
  sampled weekly and pooled monthly; one monthly pool per animal was sequenced
  while the rest are banked. Corrected 2026-08-21 from the year databases.
- Dead claim: "WF24 is genuinely 30% richer at matched depth." That difference is
  carried by sequencing run, not cohort. Corrected 2026-08-21 by
  results/rarefaction_95.tsv and logs/step7b.*.log.

NOT YET TRUSTWORTHY
- Pfam (results/pfam_hs) was still running when this block was written, jobs
  27712560 and 27712562. Quote no Pfam coverage number until 619 .done files exist.
- The composition of the unannotated 51.73% is unknown. The taxonomy result covers
  annotated proteins only and must not be extended to the dark fraction.

NEXT STEP
Union of eggNOG and Pfam coverage, then the four-way known-with-domain,
known-without-domain, genomic-unknown, environmental-unknown split, then abundance
on the support-filtered catalog.

## 2026-08-23 abundance groundwork, eggNOG+Pfam union staged

VERIFIED (file path behind each number)
- Supported set materialized: 2,069,453 reps (>=3 members, >=2 metagenomes, 95%
  tier), cross-checked exact against 40,550,595 lines and 6,182,117 clusters.
  results/supported_reps_95.txt, results/supported_reps_95.summary.txt
- WF22 CDS renamed for all 44 stems, ID sets identical to v1 input_pep (md5,
  count 0/44 mismatch). catalog/input_cds now 89/89. logs/step16e.log
- Mapping reference extracted, 2,069,453/2,069,453 found, 0 missing:
  catalog/db/LsPS_CDS_95_supported.fasta, 1,750,474,248 bp, mean 845.9
  (results/supported_cds_lengths.tsv). Index log ref_seq_len is exactly 2x this.
- Reads resolved via Fecal/read_manifest.csv: 89/89 matched, all files exist,
  89 distinct R1 and 89 distinct R2. results/reads_89.tsv
- WF24 treatment codes: 1-6 have 5 animals, code 7 has 6 (total 36).
  metadata/wf23_wf24_sequenced_metadata.tsv
- eggnog2 complete at 619/619 annotation files. 64 empty emappertmp dirs removed.
- Pfam rate measured from .done mtimes: 68.4 chunks/h, ETA ~07:00 Aug 23.

SUPERSEDES
- Dead claim (README sample table): "WF24 across 7 treatment codes with 5
  animals each." Actual: codes 1-6 x5, code 7 x6. Corrected 2026-08-23 from
  metadata/wf23_wf24_sequenced_metadata.tsv. README edit still to apply.

NOT YET TRUSTWORTHY
- Any Pfam number: quote nothing until 619 .done files exist.
- union_rep_level.tsv does not exist yet; job 27713842 runs ~08:10 with a
  619-done guard and an eggNOG cross-check (must equal 2,984,013).
- bwa-mem2 index (job 27713810, v2.3 avx2) and pilot counts (step16m, one
  sample, primary and MAPQ>=10 criteria) pending; counting criterion for the
  89-array is undecided until the pilot is reviewed.

NEXT STEP
Confirm 619 done, read 12d union table, write 12e family-tier union rollup
(cross-checks: 1,022,162 supported families, 565,769 eggNOG-dark), review
pilot, choose counting criterion, launch the 88-sample array.

## 2026-08-23 (afternoon) union verified, rep classes, GU/EU screen launched

VERIFIED (file path behind each number)
- Pfam complete: 619/619 chunks, all 8 tasks COMPLETED. results/pfam_hs
- Rep-level union (results/union_rep_level.tsv, eggNOG cross-check 2,984,013
  passed): of 6,182,117 reps, eggNOG 2,984,013 (48.27%), Pfam 2,713,156
  (43.89%, --cut_ga), union 3,227,284 (52.20%), neither 2,954,833 (47.80%).
  Supported: known 1,483,238 (71.67%), dark 586,215 (28.33%).
- Family-tier union (results/family_union_dark.tsv, anchors 2,922,537 and
  1,022,162 and 565,769 all reproduced): supported union-dark 525,557 =
  51.42%, down from 55.35% eggNOG-only. All families 72.51% from 76.01%.
- Per-rep classes (results/rep_classes.tsv, PASS): supported K 1,344,403,
  KWP 138,835, U 586,215. DUF/UPF-only K 76,589 (64,241 with eggNOG,
  12,348 Pfam-only); strict K 1,267,814.
- Mapping array 27716160 healthy; first three samples, all WF22, mapped
  fraction 0.9643 to 0.9684. results/mapping_summary/
- DB inventory: ncbi/diamond/20260128 nr_cluster_seq.dmnd 182G, dbinfo
  470,748,714 seqs, readable by diamond 2.1.24; nr.dmnd 487G Jan 2026;
  UniRef90 2022 67G readable. 2025_03 uniref100.dmnd anomalously small
  (116,964,519 seqs), excluded. kaiju group_fastas are nr-derived.

INCIDENT
- step13d glob-selected uniprot_sprot.dmnd and wrote SwissProt hits under
  the filename step13_uniref90_hits.tsv; relabeled step13_sprot2025_hits.tsv
  (or removed; count in logs/step13d.log). Byproduct, correctly labeled:
  9,633/586,215 supported-U reps hit SwissProt 2025 at E<=1e-5.
  Rule: boundary-defining databases are hard-coded, never glob-selected.

NOT YET TRUSTWORTHY
- The cross-catalog comparison (GMGC, OMRGC, AGNOSTOS) was written against
  eggNOG-only 55.35%; do not transfer it to 51.42% until criteria matched.
- step13f launched, results/step13_nrclust20260128_hits.tsv pending; no
  GU/EU number exists yet.
- Counting criterion (primary vs MAPQ>=10) open until matrix assembly.

NEXT STEP
13f to completion, EU = no hit at a stated criterion, residue verified vs
full nr.dmnd; array to 89, per-sample QC table, count matrix, contrasts.

## 2026-08-24 abundance complete, dark fraction quantified

VERIFIED
- Mapping array 89/89, no failures. Mapped fraction by cohort: WF22 0.9587
  (n=44), WF23 0.9634 (n=9), WF24 0.9671 (n=36). results/mapping_qc_89.tsv
- Count matrices built, column sums match summaries 89/89:
  results/count_matrix_primary.tsv.gz, count_matrix_mapq10.tsv.gz
  Genes with any count: primary 2,069,187; mapq10 2,042,647 of 2,069,453.
- Prevalence (results/prevalence_primary.tsv.gz): U genes in 21+ of 89
  samples 82.9%; in all 89: K 27,453, KWP 1,299, U 3,471. Only 0.7% of U
  in 5 or fewer samples.
- Dark abundance (results/dark_abundance_by_sample.tsv): U = 7.03% of raw
  mapped reads, 19.28% length-normalized. Stable across cohorts (WF22 7.11,
  WF23 6.29, WF24 7.11). Per-sample range 4.39% to 16.13%.
- Mean CDS length by class (LC_ALL=C join, all 2,069,453): K 1,045 bp,
  KWP 795, U 402.
- Completeness by class (results/completeness_by_class.tsv, flags for
  2,069,453/2,069,453): complete ORFs K 79.5%, KWP 69.7%, U 66.9%.
  Mean length of COMPLETE only: K 1,093, KWP 847, U 405.
  Conclusion: U genes are small complete proteins, not fragments.

NOT YET TRUSTWORTHY
- step13f still running (10h+); no GU/EU number yet.
- Counting criterion still primary by default; mapq10 matrix exists as
  sensitivity check, not yet compared.

NEXT STEP
Within-cohort treatment contrasts (Step 17); 13f to completion.

## 2026-08-24 dark fraction characterized, treatment design corrected

VERIFIED
- Mapping 89/89. Mapped fraction WF22 0.9587, WF23 0.9634, WF24 0.9671.
  results/mapping_qc_89.tsv
- Count matrices, column sums match summaries 89/89. Genes with any count:
  primary 2,069,187; mapq10 2,042,647 of 2,069,453.
  results/count_matrix_primary.tsv.gz, count_matrix_mapq10.tsv.gz
- Prevalence: 82.9% of U genes in 21+ of 89 samples; in all 89: K 27,453,
  KWP 1,299, U 3,471. results/prevalence_primary.tsv.gz
- Dark abundance: U = 7.03% of raw reads, 19.28% length-normalized.
  WF22 7.11, WF23 6.29, WF24 7.11. results/dark_abundance_by_sample.tsv
- Mean CDS length (LC_ALL=C join, all 2,069,453): K 1,045, KWP 795, U 402.
- Completeness (flags for 2,069,453/2,069,453): complete ORFs K 79.5%,
  KWP 69.7%, U 66.9%. Complete-only mean bp: K 1,093, KWP 847, U 405.
  So U genes are small complete proteins, not fragments.
  results/completeness_by_class.tsv, completeness_reps.tsv
- AntiFam v6.0 screen of the 586,215 supported U reps: 1,477 flagged
  (0.2520%). results/antifam_hits.tblout. Pavlopoulos 2023 found 43 of
  19,986,348 but screened a >=100-member set; criteria differ, state both.
- Cross-year core (results/core_unknown_by_years.tsv.gz): U genes in all
  3 years 497,044 (mean 419 bp, 66.2% complete, mean prevalence 45.2,
  91.18% of U reads); 2 years 69,996; 1 year 18,914 (246 bp, 1.22% of
  U reads); 0 years 261 (94 bp, no reads, assembly artifacts).
- Figures: results/fig_prevalence_by_class.png, fig_dark_abundance.png

DESIGN CORRECTION (WF22)
- WF22 is longitudinal: 44 samples from only 15 animals (14 animals x3
  months, 1 x2). metadata/wf22_design.tsv. Manuscript v1 analyzed these
  as 44 independent samples in DESeq2, so v1 treatment and egg mass
  p-values are inflated. Step 17d runs LinDA with (1|animal) plus a
  naive model to quantify the difference.
- WF22 treatment x egg mass is unbalanced: STP1710.7 occurs only in EM3
  (15 samples). Control EM1 n=2. Treatment and egg mass partly confounded.
- Aggregate class-share test (step17b) found NO treatment effect on the U
  share in any cohort (WF22 KW p=0.46/0.42, WF23 MW p=0.11/0.19,
  WF24 KW p=0.81/0.85). This tests the aggregate only, not per gene.

NOT YET TRUSTWORTHY
- step13f (nr_cluster 20260128) still running at 12h+; no GU/EU number.
- step17d LinDA results pending.
- Cross-catalog comparison text still written against eggNOG-only 55.35%.

NEXT STEP
Read 17d (mixed vs naive gene counts), finish WF23/WF24 contrasts,
13f to completion then EU call at a stated criterion.
