# Wood Frog Gut Gene Catalog: Multi-Year Extension (WF22+WF23+WF24)
# Full rewrite 2026-08-21 consolidating scaffold README + two corrections.
# APPEND ONLY from here on.

Extends LsFMGC95 (published, WF22-only, 5,055,108 genes at 95% aa id) with the
WF23 (9) and WF24 (36) sequenced fecal metagenomes. 95 WF stems total.

For all details see claude.md (state + conventions) and chat_index.md (history)
in this directory. Key facts:
- Assemblies exist for all 95 stems under
  /bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results/<stem>/
- Proteins for the 45 new stems must be predicted (Prodigal array = next step).
- Sequenced WF23/24 have no maternal or temporal axis (single library per animal;
  WF23 all EM3; WF24 33/36 EM4). They support catalog growth + treatment analyses.
  Banked longitudinal multi-EM samples exist unsequenced.
- Catalog fork (map-to-existing vs rebuild vs clusterupdate) undecided;
  map-to-existing diagnostic runs first.
- History of corrected claims is preserved in chat_index.md; nothing deleted.

## CORRECTIONS (Aug 21, 2026)
- LsFMGC95 is NOT published. It is v1 manuscript-stage, WF22-only,
  5,055,108 representative genes. Say "v1 catalog", never "published".
- Step 1a DONE: gene prediction complete for the 14 WF24 stems that lacked it
  (Prodigal v2.6.3 -p meta, same script as the v1 catalog run). All 14 succeeded
  including UHM632.41027; 10,831,320 proteins from those 14 alone.
  UHM585.41009 is a small outlier (313,220 proteins, 64 MB assembly).
- Step 1b (header renaming of all 45 to <stem>__<contig>_<n>) is the next action.
