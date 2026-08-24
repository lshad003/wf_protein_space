# Step 17h: WF23 and WF24 treatment with limma-voom. No repeated sampling in
# these cohorts, so no blocking. WF24 codes: 1 UHM520.7724, 2 UHM260.5136,
# 3 UHM516.7697, 4 Control (reference), 5 STP1717.1, 6 STP1710.7,
# 7 UHM207.4505. Codes 5 and 6 repeat the WF22 strains.
suppressMessages({library(limma); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"

meta <- read.delim(paste0(W,"/metadata/samples_89_treatment.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
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
      paste0(W,"/results/limma_", tag, "_", gsub("[^A-Za-z0-9]","_",cf), ".csv"))
  }
}

m23 <- meta[meta$cohort=="WF23", ]
m23$treatment <- relevel(factor(m23$treatment), ref="Control")
cat("WF23 groups:", table(m23$treatment), "\n")
run("WF23", m23, 8, "")

m24 <- meta[meta$cohort=="WF24", ]
m24$treatment <- relevel(factor(m24$treatment), ref="4")
m24$egg_mass <- factor(m24$egg_mass)
cat("\nWF24 groups:", levels(m24$treatment), "n =", table(m24$treatment), "\n")
cat("WF24 egg masses:", levels(m24$egg_mass), "n =", table(m24$egg_mass), "\n")
run("WF24", m24, 30, "")
cat("\nDONE step17h\n")
