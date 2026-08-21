# Wood Frog Gut Protein Space (WF22 + WF23 + WF24)

Extension of the v1 wood frog gut gene catalog (LsFMGC95) into a multi-year,
multi-tier protein-space resource. Written 2026-08-21; APPEND-ONLY below this
line except for documented full rewrites.

## Status: Step 1 complete, Step 2 in progress

## What this project is
LsFMGC95 v1 (NOT published; manuscript stage) clustered WF22 proteins at four
identity tiers. This project adds the 45 sequenced WF23 + WF24 samples and
reframes the work as a protein-space catalog: size, growth per cohort, unknown
("dark") fraction, and novel protein families.

## Questions
- Q1 How much of the wood frog gut protein space did v1 already capture?
- Q2 What is the right multi-year catalog: incremental update or rebuild?
- Q3 Does Basidiobolus treatment shift gene abundance in the new cohorts?
- Q4 How do cohorts differ, separable from batch/platform? (descriptive only)
- Q5 What are the novel/unknown proteins?
- BLOCKED: maternal and temporal replication. Not testable in sequenced
  WF23/24 (see Verified facts). Would require sequencing banked samples.

## Verified facts (all confirmed against source files)
Samples (source: ncbi-deposit/results/deposit_classification_v2.tsv,
year column is wf_year; classification is only WF-vs-wild):
- WF22 50 stems / 19 animals, longitudinal (tp-dist 1:4, 3:14, 4:1)
- WF23  9 stems /  9 animals, ONE library each (libs 40991-40999)
- WF24 36 stems / 36 animals, ONE library each (libs 41000-41034, 41909)

Design coverage of the SEQUENCED samples (45/45 joined to year databases):
- WF23 sequenced: all 9 are egg_mass 3. Treatment 5 Control vs 4 UHM520.7734.
- WF24 sequenced: 33 egg_mass 4, 3 egg_mass 2. Treatments balanced 5 x 7 codes.
- Colonies DO span egg masses (WF23: EM1 16, EM2 9, EM3 15;
  WF24: EM4 50, EM2 21, EM6 26, EM5 20, EM7 10) and WERE sampled longitudinally
  (WF23 median 38 collection dates/animal; WF24 97/107 animals >=2 weeks),
  but those samples are banked, not sequenced.

v1 catalog tiers (counted from db/LsFMGC_AA_*_rep.fasta):
- 100%: 7,545,340   95%: 5,055,108   90%: 4,547,337   50%: 2,495,918
- Build: MMseqs2 13-45111, `mmseqs cluster --min-seq-id X -c 0.8 --cov-mode 1
  --split-memory-limit 350G --kmer-per-seq 80`, on -p epyc --mem 384gb -c 96.
- Gene prediction: Prodigal v2.6.3 -p meta on scaffolds/<stem>_R.fa.gz.
- Headers: <stem>__<contig>_<n>, e.g. UHM20.35743__k141_103518_1

Step 1 outputs (this project):
- 14 WF24 stems lacked predictions; all 14 predicted successfully 2026-08-21.
- All 45 renamed to biosample-prefixed headers.
- 24,151,134 proteins across the 45 new samples
  (WF23 3,642,871 / 9 samples; WF24 20,508,263 / 36 samples).
- Per-sample counts: results/step1_protein_counts_45.tsv
- UHM585.41009 is a small outlier (313,220 proteins, 64 MB assembly).

## OPEN ISSUE 1: v1 catalog contains 2 non-wood-frog samples
The v1 build used 52 input files, not the 50 stated in the manuscript.
The extra two are UHM739.35775 and UHM740.35774, classified `wild`,
host_species `alleganiensis` (hellbender, a salamander), source `zoo`, TN,
collected 5/2023. They contributed 2,994,130 proteins (~14% of the 21.3M input).
They are absent from all WF22 metadata files and from the 44-sample analysis set.
Consequence: all catalog-descriptive numbers (input genes, representative counts,
saturation curve, annotation percentages) reflect a 52-sample build.
Abundance/DESeq2 results are structurally unaffected (only 44 WF22 samples were
mapped) but are not guaranteed identical to a clean build.
DECISION PENDING: carry forward (incremental) vs rebuild from a defined
wood-frog-only set.

## OPEN ISSUE 2: library-to-tube map missing
Nothing joins sequencing library numbers (40991+) to lab fecal biosample numbers
(WF23 7870s, WF24 8270s). Collection dates for the 45 sequenced samples are
therefore unknown. Needs the Walker group submission sheet.

## Plan
Step 0  DONE  Inventory: stems, timepoints, existing predictions, pipeline scripts
Step 1  DONE  Gene prediction for 14 missing stems + rename all 45
Step 2  NOW   Mapping-rate diagnostic: 200k-protein subsample vs LsFMGC_AA_95_rep
              GATE: decide incremental update vs rebuild (with Open Issue 1)
Step 3        Execute catalog strategy; produce multi-year catalog at 100/95/90/50
Step 4        QC layer: % complete ORFs (partial flags), AntiFam screen,
              host/diet/eukaryotic ORF classification, coverage sensitivity at -c 0.9
Step 5        Annotation + AGNOSTOS classification (eggNOG, Pfam, KEGG, CAZy,
              UniRef; DIAMOND e<=1e-5) -> K / KWP / GU / EU per cluster
Step 6        Novel family calling: >=3 members from >=2 individuals, conserved
              region >20 aa, length >=100 aa, no Pfam/eggNOG/UniRef hit,
              AntiFam-clean, dN/dS<0.5 where computable
Step 7        Structural annotation of novel reps: ESMFold -> pLDDT/pTM>0.7 ->
              Foldseek vs PDB + AFDB + ESM Atlas (expect ~14-16% yield)
Step 8        Abundance: bwa-mem2 + featureCounts, all 95 samples vs final catalog
              -> prevalence, core/accessory, cohort sharing; treatment contrasts
Step 9        Batch table for 95 stems (platform, run, depth). Parallel; required
              before any cross-year claim.
Step 10       Write-up. Figures: rarefaction per tier per cohort; cluster-size
              distribution; known/unknown stacked bar; presence/absence heatmap;
              cohort-sharing chord; novel-family showcase.

## Key paths
assemblies   /bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results/<stem>/
v1 catalog   /bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1/db/
v1 pipeline  /bigdata/stajichlab/shared/projects/Herptile/Metagenome/LsFMGC/v1/pipeline/
predictions  /bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly/
year dbs     /bigdata/stajichlab/lshad003/ncbi-deposit/wf_databases/
this project /bigdata/stajichlab/lshad003/wf_protein_space/
python       /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
cluster host lshad003@cluster.hpcc.ucr.edu

## Conventions
Scripts as heredoc blocks in scripts/, run immediately. Full paths always.
No find, no wget, no em-dashes. Python at the full conda path, never source activate.
Verification-first: run a check, paste output, then advise. Never state a count
as fact without a verified output file.
MMseqs2 MUST be pinned to 13-45111 (bare `module load mmseqs2` gives 17).
MMseqs2 jobs MUST use -p epyc (crashes with Illegal instruction on some nodes).
/bigdata was 98% full on 2026-08-21; check space before large jobs.
md files are append-only; wrong analyses are marked corrected with a date, never deleted.

## INSTRUCTIONS FOR CLAUDE CODE
Execute ONE step per session, the one the user names. After the step's checks
pass, append a dated entry to logs/step_log.md and STOP. Do not begin the next
step. If a check fails, report and STOP; do not "fix" it by changing the plan.
Never modify anything under /bigdata/stajichlab/shared/ (read-only for us).
Write all outputs under /bigdata/stajichlab/lshad003/wf_protein_space/.
Heavy work goes in an sbatch script (-p epyc for MMseqs2), never on the login node.
Do not commit data files; .gitignore covers fasta/faa/gff/bam and catalog/.
Before claiming any number, show the command that produced it.
