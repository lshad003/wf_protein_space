#!/bin/bash
# Status check: jobs since 2026-10-22 and the outputs we are waiting on.
cd /bigdata/stajichlab/lshad003/wf_protein_space
echo "== now: $(date)"
echo "== tmux sessions"
tmux ls 2>&1
echo "== still queued or running"
squeue -u lshad003 -o "%.14i %.32j %.10P %.9T %.10M %.12l %R"
echo "== jobs since 2026-10-22, count by name and state"
sacct -u lshad003 -S 2026-10-22 -X -n -P --format=JobName,State | sort | uniq -c
echo "== jobs since 2026-10-22 that did not complete"
sacct -u lshad003 -S 2026-10-22 -X -n -P --format=JobID,JobName,State,ExitCode,Elapsed,End | grep -v -e "|COMPLETED|" -e "|RUNNING|" -e "|PENDING|"
echo "== srun session 30323149 and near-miss job 30338647"
sacct -j 30323149,30338647 -X -n -P --format=JobID,JobName,State,ExitCode,Elapsed,End
echo "== hmmsearch array 30338473, task states"
sacct -j 30338473 -X -n -P --format=State,ExitCode | sort | uniq -c
echo "== newest 20 logs"
ls -lt logs | head -21
for pat in "logs/step36b*" "logs/*sharedruns*" "logs/step35*"; do
  for f in $(ls -t $pat 2>/dev/null | head -2); do
    echo "-- tail of $f"
    tail -20 "$f"
  done
done
echo "== newest files in results/step36b and results/step35"
ls -lt results/step36b 2>&1 | head -15
ls -lt results/step35 2>&1 | head -10
echo "== month shared-run outputs"
ls -l results/check_wf22_month_run.txt 2>&1
for f in results/*/month_shared_runs.tsv results/*/month1_excess_by_animal.tsv; do
  if [ -e "$f" ]; then echo "-- $f"; cat "$f"; fi
done
echo "== git"
git status --short | head -40
git log --oneline -6
