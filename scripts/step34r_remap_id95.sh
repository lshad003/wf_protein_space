#!/bin/bash
#SBATCH -J step34r_id95
#SBATCH -p stajichlab -c 16 --mem 64G --time 12:00:00
#SBATCH --array=1-63,65-89%8
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34r_remap_id95.%A_%a.log
# Step 34r: identity-filtered recount. As scripts/step16q_map_array.sh (same reference
# catalog/db/LsPS_CDS_95_supported.fasta, bwa-mem2 2.3 mem -t 16, samtools view -F 0x900, same rows of
# results/reads_89.tsv), but on fastp-trimmed reads instead of the raw FASTQs.
# Array: the 88-sample set; task 64 (UHM586.41010, excluded sample) is not run. step16q skipped task 11
# (UHM17.35739) only because outputs existed; it is mapped here.
# fastp 0.23.2 with the recorded flags (results/methods_facts.tsv:12: -w 8 --detect_adapter_for_pe --trim_poly_g
# --overrepresentation_analysis --unpaired1/--unpaired2, other settings default); -j/-h only set the report paths.
# Trimmed reads live on node-local scratch and are deleted at task end, also on failure (trap on EXIT).
# Count columns per gene: gene, primary mapped (step16q col2), MAPQ >= 10 (step16q col3), c_id95, c_id90, c_id95_q1.
#   identity = 1 - NM / (M + I + D) (M includes = and X); end-consistent clipping: a soft (or hard) clip > 3 bases is
#   allowed only at the left if the alignment starts within 3 bases of the reference start (POS <= 4), and at the right
#   if it ends within 3 bases of the reference end (reflen - end <= 3); reference lengths from
#   results/step34r/LsPS_CDS_95_supported.fasta.fai. c_id95_q1 = c_id95 and MAPQ >= 1.
# Outputs: results/step34r/counts/<S>.counts.tsv.gz, results/step34r/summary/<S>.summary.tsv (with new/old col2
#   ratio against results/counts/<S>.counts.tsv.gz), results/step34r/fastp/<S>.{json,log}.
# Exits with the main commands' status (CLAUDE.md rule).
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
REF=$WPS/catalog/db/LsPS_CDS_95_supported.fasta
FAI=$WPS/results/step34r/LsPS_CDS_95_supported.fasta.fai
O=$WPS/results/step34r
[ -s $REF.bwt.2bit.64 ] || { echo "ABORT: bwa-mem2 index missing"; exit 1; }
[ -s $FAI ] || { echo "ABORT: $FAI missing"; exit 1; }
mkdir -p $O/counts $O/summary $O/fastp
LINE=$((SLURM_ARRAY_TASK_ID+1))
IFS=$'\t' read -r S R1 R2 <<< "$(sed -n "${LINE}p" $WPS/results/reads_89.tsv)"
[ -n "$S" ] || { echo "ABORT: no row $LINE"; exit 1; }
if [ "$S" = "UHM586.41010" ]; then echo "ABORT: UHM586.41010 is excluded"; exit 1; fi
if [ -e $O/counts/$S.counts.tsv.gz ] || [ -e $O/summary/$S.summary.tsv ]; then echo "SKIP $S: outputs exist (not overwritten)"; exit 0; fi
[ -s "$R1" ] && [ -s "$R2" ] || { echo "ABORT: missing reads for $S"; exit 1; }
module load workspace/scratch
[ -n "$SCRATCH" ] || { echo "ABORT: SCRATCH not set"; exit 1; }
T=$SCRATCH/step34r_$S; mkdir -p $T
trap 'rm -rf "$T"' EXIT
echo "sample $S task $SLURM_ARRAY_TASK_ID host $(hostname) $(date)"
module load fastp/0.23.2
fastp -w 8 --detect_adapter_for_pe -i "$R1" -I "$R2" -o $T/R1.fq.gz --out2 $T/R2.fq.gz --trim_poly_g \
  --unpaired1 $T/unpair1.fq.gz --unpaired2 $T/unpair2.fq.gz --overrepresentation_analysis \
  -j $T/fastp.json -h $T/fastp.html 2> $T/fastp.log
rc=$?
cp $T/fastp.json $O/fastp/$S.json; cp $T/fastp.log $O/fastp/$S.log
if [ $rc -ne 0 ]; then echo "ABORT: fastp exit $rc"; exit 1; fi
grep -E "reads passed filter|Duplication|Insert size" $T/fastp.log
echo "fastp done $(date)"
module load bwa-mem2/2.3 samtools/1.22.1
bwa-mem2 mem -t 16 $REF $T/R1.fq.gz $T/R2.fq.gz 2> $T/bwa.err \
| samtools view -F 0x900 - \
| awk -v CF="$T/counts.tsv" -v S="$S" '
    FNR == NR { L[$1] = $2; next }
    {
      tot++
      if (int($2 / 4) % 2 == 1) { unm++; next }
      g = $3; a[g]++; if ($5 >= 10) q[g]++
      nm = -1
      for (i = 12; i <= NF; i++) if (substr($i, 1, 5) == "NM:i:") { nm = substr($i, 6) + 0; break }
      if (nm < 0) { nonm++; next }
      if (!(g in L)) { nolen++; next }
      cig = $6; m = 0; ins = 0; del = 0; refc = 0; lc = 0; rcl = 0; first = 1
      while (cig != "") {
        match(cig, /^[0-9]+/); n = substr(cig, 1, RLENGTH) + 0; op = substr(cig, RLENGTH + 1, 1); cig = substr(cig, RLENGTH + 2)
        if (op == "M" || op == "=" || op == "X") { m += n; refc += n; first = 0 }
        else if (op == "I") { ins += n; first = 0 }
        else if (op == "D") { del += n; refc += n; first = 0 }
        else if (op == "N") { refc += n; first = 0 }
        else if (op == "S" || op == "H") { if (first) lc += n; else rcl += n }
      }
      aln = m + ins + del; if (aln == 0) next
      id = 1 - nm / aln; st = $4; en = st + refc - 1
      if (lc > 3 && st > 4) next
      if (rcl > 3 && L[g] - en > 3) next
      if (id >= 0.95) { i95[g]++; if ($5 >= 1) i95q[g]++ }
      if (id >= 0.90) i90[g]++
    }
    END {
      for (g in a) printf "%s\t%d\t%d\t%d\t%d\t%d\n", g, a[g], q[g], i95[g], i90[g], i95q[g] > CF
      printf "%s\t%d\t%d\t%.4f\t%d\t%d\n", S, tot, tot - unm, (tot - unm) / tot, nonm, nolen
    }' $FAI - > $T/summary.tsv
ps=("${PIPESTATUS[@]}"); rc=0
for x in "${ps[@]}"; do if [ "$x" -ne 0 ]; then rc=$x; break; fi; done
echo "pipeline exit codes bwa-mem2/samtools/awk: ${ps[*]}"
if [ $rc -ne 0 ]; then tail -5 $T/bwa.err; echo "ABORT: mapping pipeline exit $rc"; exit 1; fi
NEW=$(awk '{s += $2} END {print s + 0}' $T/counts.tsv)
OLD=$(zcat $WPS/results/counts/$S.counts.tsv.gz | awk '{s += $2} END {print s + 0}')
RATIO=$(awk -v n=$NEW -v o=$OLD 'BEGIN {if (o > 0) printf "%.4f", n / o; else print "NA"}')
gzip -c $T/counts.tsv > $O/counts/$S.counts.tsv.gz
rc=$?
if [ $rc -ne 0 ]; then echo "ABORT: gzip exit $rc"; rm -f $O/counts/$S.counts.tsv.gz; exit 1; fi
awk -v n=$NEW -v o=$OLD -v r=$RATIO 'BEGIN {OFS = "\t"} {print $0, o, n, r}' $T/summary.tsv > $O/summary/$S.summary.tsv
rc=$?
echo "columns: sample total_records mapped frac_mapped mapped_without_NM mapped_without_ref_length old_col2 new_col2 new_over_old"
cat $O/summary/$S.summary.tsv
echo "DONE $S $(date)"
exit $rc
