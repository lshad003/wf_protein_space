# claude.md - WF Protein Space (multi-year catalog extension)
# READ THIS FIRST in every new chat for this project.
# Last full rewrite: 2026-08-21. After this date, updates are APPEND-ONLY (cat >>).

## WHO / WHAT
Leila Shadmani (lshad003), PhD student, Stajich Lab, UC Riverside.
Project: extend the published wood frog gut gene catalog LsFMGC95 (WF22 only)
with WF23 and WF24 collections.

## WORKING CONVENTIONS (never violate)
- Scripts as ready-to-run heredoc blocks, saved under
  /bigdata/stajichlab/lshad003/wf_protein_space/scripts/ and run immediately.
- Python at full conda path: /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
  Never source activate.
- No find, no wget, no em-dashes, no set -euo pipefail. Use ls output pasted back.
- Verification-first: Leila runs checks on the cluster and pastes output BEFORE
  Claude advises. Never state counts as fact without a verified output file.
- Always write FULL paths in every command. Never ask Leila to edit manually.
- md files: cat >> append only. cat > only if the block contains the complete
  prior content. Wrong analyses get marked corrected with date, never deleted.
- One topic per chat. Title format: [CATEGORY] topic - detail
  (ANALYSIS, FIGURE, AUDIT, METHODS, RESULTS, STATS, CODE, SETUP).
- End of every significant chat: append entry to chat_index.md, append any new
  facts here, remind Leila to re-upload both files to the Claude project sidebar.

## LsFMGC95 (published v1, DO NOT change these numbers)
- WF22: 50 fecal shotgun metagenomes, Jan 2022 egg masses, one pond White Co TN,
  EM1/EM2/EM3, lab-reared, Basidiobolus (STP1717.1, STP1710.7) vs control,
  3 monthly timepoints. 44/50 post-QC (EM1=8, EM2=9, EM3=27).
- Pipeline: fastp; metashot/mag-illumina v2.2.0 (MEGAHIT, contigs >=1500bp);
  Prodigal -p meta; MMseqs2 v13, 95% aa id, -c 0.8 --cov-mode 1 --kmer-per-seq 80
  = 5,055,108 representative genes. Headers: UHM20.35743__k141_103518_1 style.
- Main result: maternal origin (egg mass) dominates gene abundance and
  strengthens over development.

## VERIFIED STATE OF THE EXTENSION (as of 2026-08-21)
Sample lists (source: /bigdata/stajichlab/lshad003/ncbi-deposit/results/deposit_classification_v2.tsv,
year column is wf_year, NEVER key on classification which is only WF-vs-wild):
- WF22 50 stems / 19 animals (longitudinal: tp-dist 1:4, 3:14, 4:1)
- WF23  9 stems /  9 animals, one library each (libs 40991-40999)
- WF24 36 stems / 36 animals, one library each (libs 41000-41034 + 41909)
- Stem lists: /bigdata/stajichlab/lshad003/wf_protein_space/metadata/wf22_stems.txt (50),
  wf23_stems.txt (9), wf24_stems.txt (36)

Sequenced-sample metadata (verified join, 45/45 matched):
- /bigdata/stajichlab/lshad003/wf_protein_space/metadata/wf23_wf24_sequenced_metadata.tsv
  columns: stem, animal, wf_year, egg_mass, treatment, tank, metamorphosis_date
- WF23 sequenced: ALL 9 animals egg_mass 3. Treatment 5 Control vs 4 UHM520.7734.
- WF24 sequenced: 33/36 egg_mass 4, 3 egg_mass 2 (UHM590, UHM610, UHM624).
  Treatments balanced 5 per code x 7 codes.
- CONSEQUENCE: sequenced WF23/24 = catalog growth + treatment contrasts ONLY.
  No maternal validation, no developmental axis, without sequencing banked samples.

Colonies were sampled longitudinally (banked in freezer, NOT sequenced):
- WF23: 40 animals, median 38 fecal collection dates each, Mar-Jul 2023,
  egg masses 1/2/3, Control + 4 Basidiobolus isolates.
- WF24: 144 animals, weekly pooled feces, 97/107 with >=2 weeks,
  egg masses 2/4/5/6/7, 7 treatment codes.
- Year databases: /bigdata/stajichlab/lshad003/ncbi-deposit/wf_databases/ (3 xlsx)

Infrastructure:
- Assemblies EXIST for all 95 WF stems:
  /bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results/<stem>/
  (megahit/, scaffolds/, bins/, metabat2/, qc/, unbinned/)
- Per-sample predicted PROTEINS DO NOT EXIST for the 45 new stems.
  NEXT STEP: Prodigal -p meta SLURM array over 45 stems, headers biosample-prefixed
  <stem>__<contig>_<n> to match LsFMGC95 scheme.

## OPEN DECISIONS AND BLOCKERS
1. Catalog strategy fork: map-to-existing vs rebuild vs MMseqs2 clusterupdate.
   FIRST run map-to-existing diagnostic: WF23/24 mapping rate against LsFMGC95.
   Rebuild renumbers every gene ID and invalidates all verified WF22 outputs.
2. BLOCKER library-to-tube map: nothing joins sequencing library numbers (40991+)
   to lab fecal biosample numbers (WF23 7870s, WF24 8270s). Collection dates for
   the 45 sequenced samples unknown. Ask Walker group for the submission sheet.
3. Confounders: year confounded with platform/batch/depth (WF22 spanned
   NovaSeq + AVITI). Need per-sample table: year, platform, run, depth.
4. Possible future sequencing of banked samples (other egg masses, other weeks)
   would restore maternal + temporal axes. Needs Jason discussion.

## MANUSCRIPT LANGUAGE RULES (inherited from WF22 audits, still binding)
- froglets not tadpoles for post-metamorphic animals.
- "functional" language ONLY for KEGG/Pfam/CAZy-annotated analyses;
  DESeq2 counts, distances, correlations = "gene abundance" only.
- Pfam is annotation not clustering; 95% clusters are the analytical units.

## CLUSTER LOGIN (appended 2026-08-21)
scp/ssh host: lshad003@cluster.hpcc.ucr.edu

## CORRECTION (Aug 21, 2026): LsFMGC95 status
Earlier text called LsFMGC95 "published". It is NOT published: it is v1
manuscript-stage, WF22-only, 5,055,108 representative genes, drafted and
circulated but not submitted/accepted.
Implication: no published numbers are locked by a rebuild. Incremental
clustering is still preferred, but for a different reason: it preserves gene
IDs already used in v1 figures, DESeq2 outputs, and the per-sample annotation
tables shared with Joshua Phillips, avoiding redoing verified work.
Language rule going forward: say "v1 manuscript" or "v1 catalog", never
"published catalog".
