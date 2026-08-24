# CLAUDE.md - operating rules for this project
# State lives in CHATINDEX.md, the analysis record in README.md, history in
# PROJECT_LOG.md. Read CHATINDEX.md first; do not re-derive settled questions.

## Cluster
- python: /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3, never source activate
- partitions: epyc and stajichlab are 30 days, short is 2 hours only
- MMseqs2 pinned to 13-45111; bare `module load mmseqs2` gives 17
- MMseqs2 runs only on -p epyc; other nodes fail with an illegal instruction
- large reference databases must be staged to node-local /scratch or concurrent
  tasks stall; this applies to eggNOG (48 GB) and Pfam (4.7 GB)
- emapper --annotate_hits_table returns zero rows in 2.1.9 and 2.1.7; every
  chunk must run end to end
- Pfam: use hmmsearch, not hmmscan. hmmscan measured 78 min per 10k chunk
  against about 6 min for hmmsearch. Column order differs between them.

## Method rules
- coverage mode is stated explicitly for every search; --cov-mode 2 for searches
- cluster-level claims are made on support-filtered clusters
- cohort comparisons are made at matched sampling effort and checked against
  sequencing run
- no count is stated without the output file behind it
- a count that depends on a threshold is reported with its criterion

## Working conventions
- scripts as ready-to-run heredocs into scripts/, full paths always
- no find, no wget, no set -euo pipefail, no em-dashes
- .md files are append-only except on an explicit request for a rewrite
- data files are never committed; .gitignore is a whitelist
- LinDA (MicrobiomeStat) allocates an n_features^2 matrix and fails at this
  scale: 368,236 genes asked for 1010.3 Gb. Use limma-voom instead, with
  duplicateCorrelation for repeated sampling.
- WF22 is longitudinal: 44 samples from 15 animals across 3 months. Any
  WF22 model must block on animal or the p-values are inflated.
