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
