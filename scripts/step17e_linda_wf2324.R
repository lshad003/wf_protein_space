# Step 17e: LinDA treatment models for WF23 and WF24.
# No repeated sampling in these cohorts, so no random effect.
suppressMessages(library(MicrobiomeStat))
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"

meta <- read.delim(paste0(W,"/metadata/samples_89_treatment.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
cat("reading count matrix...\n")
cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)

run <- function(tag, mt, formula, minprev) {
  cn <- cnt[, mt$stem, drop=FALSE]
  keep <- rowSums(cn > 0) >= minprev
  cn <- cn[keep, , drop=FALSE]
  cat("\n=== ", tag, " ===\n", sep="")
  cat("samples:", ncol(cn), " genes after prevalence filter:", nrow(cn), "\n")
  r <- linda(feature.dat=cn, meta.dat=mt, formula=formula,
             feature.dat.type="count", prev.filter=0, mean.abund.filter=0,
             is.winsor=FALSE, n.cores=4, verbose=FALSE)
  for (v in names(r$output)) {
    o <- r$output[[v]]
    cat(sprintf("  %-24s padj<0.05: %6d   +|LFC|>1: %6d\n", v,
      sum(o$padj < 0.05, na.rm=TRUE),
      sum(o$padj < 0.05 & abs(o$log2FoldChange) > 1, na.rm=TRUE)))
    write.csv(o[order(o$padj), ],
      paste0(W,"/results/linda_", tag, "_", gsub("[^A-Za-z0-9]","_",v), ".csv"))
  }
}

m23 <- meta[meta$cohort=="WF23", ]
m23$treatment <- relevel(factor(m23$treatment), ref="Control")
run("WF23", m23, "~ treatment", 9)

m24 <- meta[meta$cohort=="WF24", ]
m24$treatment <- relevel(factor(m24$treatment), ref="4")
cat("\nWF24 treatment levels:", levels(m24$treatment), "\n")
run("WF24", m24, "~ treatment", 30)
cat("\nDONE step17e\n")
