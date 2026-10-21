# Step 34s: one side-by-side table of the read-count numbers, old matrix vs identity-filtered matrices.
# Old = results/step34s/old/ (byte-identical to the original outputs, results/step34s/repro_old_check.tsv);
#   c_id95 (primary), c_id90 and c_id95_q1 (sensitivity) = results/step34s/<label>/ (scripts/step34s_run.sh).
#   Sensitivity labels have PERMANOVA and limma only (part rsens); other cells are NA.
# Rows: read share, prevalence and three-cohort presence (step34b2_catalog_stats_final.tsv, strict_final rows);
#   read retention (step34k_summary.tsv Part D); Step 34e PERMANOVA and limma (step34e_<set>_permanova/limma.tsv);
#   Step 34l months 2-3 PERMANOVA (step34l_m23_permanova.tsv); Fig. 5a,b PCoA axis % (fig_pcoa_wf22.tsv).
# Output: results/step34s/side_by_side.tsv; refuses to overwrite.
import os, csv
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; S = W + "/results/step34s/"
OUT = S + "side_by_side.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
LABS = ["old", "c_id95", "c_id90", "c_id95_q1"]
def rd(lab, f):
    p = S + lab + "/" + f
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else None
T = {}  # key -> {lab: value}; key = (section, item, file)
def put(sec, item, f, lab, v): T.setdefault((sec, item, f), {})[lab] = v
for lab in LABS:
    st = rd(lab, "step34b2_catalog_stats_final.tsv")
    if st:
        for r in st:
            if r["set"] != "strict_final": continue
            for col, nm in (("pct_reads_raw_pooled", "% of reads (pooled)"), ("pct_reads_lennorm_pooled", "% of reads, length-normalised (pooled)"),
                            ("mean_prev_ge1read", "mean prevalence, >= 1 read (of 88)"), ("mean_prev_ge5reads", "mean prevalence, >= 5 reads (of 88)"),
                            ("pct_3cohorts_ge1read", "% in all three cohorts, >= 1 read"), ("pct_3cohorts_ge5reads", "% in all three cohorts, >= 5 reads")):
                put("read share / prevalence (strict_final)", "%s: %s" % (r["class"], nm), "step34b2_catalog_stats_final.tsv", lab, r[col])
    sm = rd(lab, "step34k_summary.tsv")
    if sm:
        for r in sm:
            if r["section"] == "D" and ("retention" in r["metric"]):
                put("read retention (Fig. 1d)", r["metric"], "step34k_summary.tsv", lab, r["value"])
    for s in ("strict_final", "strict_final_EU"):
        pm = rd(lab, "step34e_%s_permanova.tsv" % s)
        if pm:
            for r in pm:
                if r["note"].startswith("untestable"): continue
                it = "%s | %s | %s | %s | %s" % (s, r["cohort"], r["metric"], r["test"], r["term"])
                v = ("F %s p %s" % (r["F"], r["p"])) if r["R2"] == "NA" else ("R2 %s p %s" % (r["R2"], r["p"]))
                put("PERMANOVA (Step 34e)", it, "step34e_%s_permanova.tsv" % s, lab, v)
        lm = rd(lab, "step34e_%s_limma.tsv" % s)
        if lm:
            for r in lm:
                it = "%s | %s | %s | %s" % (s, r["cohort"], r["model"], r["coef"])
                v = r["padj_lt_0.05"] if r["coef"] == "within_animal_correlation" else "%s of %s" % (r["padj_lt_0.05"], r["genes_tested"])
                put("limma-voom (Step 34e; padj < 0.05 of genes tested; consensus correlation)", it, "step34e_%s_limma.tsv" % s, lab, v)
    m23 = rd(lab, "step34l_m23_permanova.tsv")
    if m23:
        for r in m23:
            put("PERMANOVA months 2-3 (Step 34l)", "%s | %s | %s | %s" % (r["set"], r["metric"], r["test"], r["term"]), "step34l_m23_permanova.tsv", lab,
                "R2 %.5f p %s" % (float(r["R2"]), r["p"]))
    pc = rd(lab, "fig_pcoa_wf22.tsv")
    if pc:
        for s in ("strict_final", "strict_final_EU"):
            r = next(x for x in pc if x["set"] == s)
            put("PCoA (Fig. 5a,b)", "%s: PCo1 / PCo2 %%" % s, "fig_pcoa_wf22.tsv", lab, "%.1f / %.1f" % (float(r["pct1"]), float(r["pct2"])))
with open(OUT, "w") as o:
    o.write("section\titem\t" + "\t".join(LABS) + "\tfile (in results/step34s/<label>/)\n")
    for (sec, it, f), d in T.items():
        o.write("\t".join([sec, it] + [d.get(l, "NA") for l in LABS] + [f]) + "\n")
print("rows %d; wrote %s" % (len(T), OUT))
