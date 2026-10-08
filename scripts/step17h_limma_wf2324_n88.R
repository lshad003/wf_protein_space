# Step 17h n88: WF23 part of step17h_limma_wf2324.R on the 88-metagenome
# matrix and supported set. The WF24 part is replaced by step17i_limma_wf24_n34_n88.R.
# No repeated sampling in WF23, so no blocking.
suppressMessages({library(limma); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"

meta <- read.delim(paste0(W,"/metadata/samples_88_treatment.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary_n88.tsv.gz")),
                  row.names=1, check.names=FALSE)

run <- function(tag, mt, minprev, extra) {
  cn <- cnt[, mt$stem, drop=FALSE]
  keep <- rowSums(cn > 0) >= minprev
  cn <- cn[keep, , drop=FALSE]
  cat("\n=== ", tag, " ===\n", sep="")
  cat("samples:", ncol(cn), " genes present in >=", minprev, ":", nrow(cn), "\n")
  d <- calcNormFactors(DGEList(counts=as.matrix(cn)), method="TMM")
  des <- model.matrix(as.formula(paste("~ treatment", extra)), data=mt)
  fit <- eBayes(lmFit(voom(d, des), des))
  for (cf in colnames(coef(fit))) {
    if (cf == "(Intercept)") next
    t <- topTable(fit, coef=cf, number=Inf, sort.by="none")
    cat(sprintf("  %-24s padj<0.05: %6d   +|LFC|>1: %6d\n", cf,
      sum(t$adj.P.Val < 0.05, na.rm=TRUE),
      sum(t$adj.P.Val < 0.05 & abs(t$logFC) > 1, na.rm=TRUE)))
    write.csv(t[order(t$adj.P.Val), ],
      paste0(W,"/results/limma_", tag, "_", gsub("[^A-Za-z0-9]","_",cf), "_n88.csv"))
  }
}

m23 <- meta[meta$cohort=="WF23", ]
m23$treatment <- relevel(factor(m23$treatment), ref="Control")
cat("WF23 groups:", table(m23$treatment), "\n")
run("WF23", m23, 8, "")
cat("\nDONE step17h_n88\n")
