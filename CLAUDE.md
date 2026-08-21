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
- WF23  9 stems /  9 animals, one monthly-pool metagenome per animal (timepoint IDs 40991-40999)
- WF24 36 stems / 36 animals, one monthly-pool metagenome per animal (timepoint IDs 41000-41034, 41909)

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

## WF22 analysis set (v1 manuscript)
50 WF22 samples were collected and used to build the v1 catalog; 44 passed
review and carried into abundance and statistical analyses.

Six exclusions:

| sample | reason |
|--------|--------|
| UHM56.10839  | anomalously high Basidiobolus reads in a control (3.09%) |
| UHM102.35765 | post-sequencing quality review |
| UHM27.10829  | animal died before completing the 3-month experiment |
| UHM33.10831  | animal died before completing the 3-month experiment |
| UHM43.10836  | animal died before completing the 3-month experiment |
| UHM44.10837  | animal died before completing the 3-month experiment |

Final analysis set: 44 samples (EM1 8, EM2 9, EM3 27). EM3 is larger because
it received both Basidiobolus strains (STP1717.1 and STP1710.7) plus controls.

Note for the extension: catalog construction and abundance analysis use
different sample sets. Any catalog-level count refers to the build set;
any abundance or statistical result refers to the analysis set.


## Longitudinal sampling in WF23 and WF24 (clarification)
The sequenced animals WERE sampled repeatedly over time. Verified from the year
databases: each of the 9 sequenced WF23 animals has 15 to 42 dated fecal
collections with lab biosample IDs; each of the 36 sequenced WF24 animals has
5 to 7, held in a sheet named `pooled_samples`.

WF22 practice was to pool fecal material by month before sequencing (see the
`pool_month` tables in the v1 directory), and the WF24 sheet name indicates the
same practice, so each WF23/WF24 library most likely represents pooled material
spanning several collection dates rather than a single moment.

What this means analytically: WF22 has 3 to 4 libraries per animal and therefore
supports within-animal comparison across development; WF23 and WF24 have one
library per animal, so they do not, regardless of how many collections went into
that library.

Not yet verifiable: which specific lab collections were pooled into each
sequencing library. Requires the sequencing submission sheet that maps library
numbers (40991+) to lab biosample numbers (WF23 8070s-18000s, WF24 8270s).


## Sampling and pooling design (WF22 template)
Fecal samples were collected weekly and pooled by month; each sequenced
metagenome is one month's pool for one animal. WF22 therefore has 3 metagenomes
per animal across the 3-month experiment (e.g. UHM20 -> 10828, 35743, 35744),
which is what supports the within-animal developmental comparison.

WF23 and WF24 have one metagenome per animal, so one month's pool per animal was
sequenced while the remaining monthly pools stay banked. Which month each
sequenced pool represents is not yet known; it requires the submission sheet
mapping timepoint IDs (40991+) to lab biosample IDs. This matters because if
sequenced pools come from different months across animals, month is an
uncontrolled variable in the WF24 treatment comparison.

### Step 2 - mapping-rate diagnostic (2026-08-21)
200k-protein random subsample of the 45 new samples searched against v1 tiers,
MMseqs2 13-45111, -s 4, e<=1e-5, -c 0.8.

A WF22 positive control (proteins that built the catalog, so expected ~100%)
showed that `--cov-mode 1` and `--cov-mode 0` systematically miss short proteins
(control only 61.4% and 60.9%). `--cov-mode 2` (coverage of the query) recovers
the control at 98.25% and is the correct setting for this comparison.

Verified rates under cov-mode 2:

| search | rate |
|--------|------|
| control: WF22 vs 95% tier | 98.25% |
| new 45 vs 95% tier (gene level) | 46.15% |
| new 45 vs 50% tier (family level) | 66.66% |

Interpretation: roughly 54% of WF23/WF24 proteins are novel at gene level and
33% at family level. Gene-level space is far from saturated while family-level
space is better covered, matching the UHGP pattern.
Scripts: scripts/step2c_diagnostic_epyc.sh, step2d_control_and_tiers.sh,
step2e_covmode2.sh
Note: MMseqs2 crashes with Illegal instruction on non-epyc nodes; always use -p epyc.

## Sample structure, verified
Fecal samples were collected weekly and pooled by month. Each sequenced
metagenome is one monthly pool for one animal; the stem is
<animal>.<timepoint-ID>, e.g. UHM20.10828, UHM20.35743, UHM20.35744 are three
consecutive monthly pools for animal UHM20.

- WF22: 50 metagenomes from 19 animals, multiple monthly pools per animal
  (4 animals with 1, 14 with 3, 1 with 4). Supports within-animal comparison
  across development.
- WF23: 9 metagenomes from 9 distinct animals, one monthly pool each,
  all from the same month.
- WF24: 36 metagenomes from 36 distinct animals, one monthly pool each.

Remaining monthly pools for WF23 and WF24 animals are banked, not sequenced.
WF24 metamorphosis dates span 2024-04-07 to 2024-05-19, so time since
metamorphosis at collection varies across animals.
