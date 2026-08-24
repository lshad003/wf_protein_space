#!/bin/bash
# Which file holds ORF completeness, and how step9 defined it.
BASE=/bigdata/stajichlab/lshad003/wf_protein_space
echo "=== step9 script: inputs and outputs ==="
grep -nE 'partial|complete|open\(|results/' $BASE/scripts/step9*.* 2>/dev/null | head -20
echo ""
echo "=== candidate completeness files in results ==="
ls -l $BASE/results | grep -iE 'complet|partial|orf'
