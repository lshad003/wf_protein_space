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
