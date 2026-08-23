# CLAUDE.md - operating rules for this project

State lives in CHATINDEX.md (append-only, newest block wins, a claim named in a
SUPERSEDES line is dead). The analysis record is README.md. History is
PROJECT_LOG.md. Read CHATINDEX.md first and do not re-derive settled questions.

## Cluster facts
- python: /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3, never source activate
- partitions: epyc and stajichlab allow 30 days, short caps at 2 hours
- MMseqs2 pinned to 13-45111; bare `module load mmseqs2` gives version 17
- MMseqs2 runs only on -p epyc; other nodes fail with an illegal instruction
- large reference databases must be staged to node-local /scratch or concurrent
  tasks stall: eggNOG is 48 GB, Pfam 38.2 is 4.7 GB
- emapper --annotate_hits_table returns zero rows in both 2.1.9 and 2.1.7, so
  every chunk runs end to end
- Pfam uses hmmsearch, not hmmscan: 78 minutes per 10,000-protein chunk against
  about 6 minutes. Their table column orders differ.

## Method rules
- coverage mode is stated explicitly for every search; searches use --cov-mode 2
- cluster-level claims are made on support-filtered clusters, not raw clusters
- cohort comparisons are made at matched sampling effort and checked against
  sequencing run before interpretation
- no number is stated without the current output file behind it
- any count that depends on a threshold is reported with its criterion
- a selection criterion is never reported as a finding

## Working conventions
- scripts as ready-to-run heredocs written into scripts/, full paths always
- no find, no wget, no set -euo pipefail, no em-dashes
- .md files are append-only unless a full rewrite is explicitly requested
- data files are never committed; .gitignore is a whitelist
- one topic per chat; end each substantive chat with an append block for
  CHATINDEX.md and a one-line append to PROJECT_LOG.md
