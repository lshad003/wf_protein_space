# Step 17g: WF22 treatment effect with limma-voom.
# Repeated sampling of the same animal handled with duplicateCorrelation,
# the standard approach for a blocked design. Memory is linear in genes.
suppressMessages({library(limma); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
MINPREV <- 35

meta <- read.delim(paste0(W,"/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
meta$treatment <- relevel(factor(meta$treatment), ref="control")
meta$egg_mass <- factor(meta$egg_mass)
meta$month <- factor(meta$month)
meta$animal <- factor(meta$animal)
cat("samples:", nrow(meta), " animals:", nlevels(meta$animal), "\n")

cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)
cnt <- cnt[, meta$stem, drop=FALSE]
keep <- rowSums(cnt > 0) >= MINPREV
cnt <- cnt[keep, , drop=FALSE]
cat("genes present in >=", MINPREV, "of 44:", nrow(cnt), "\n")

d <- DGEList(counts=as.matrix(cnt))
d <- calcNormFactors(d, method="TMM")

report <- function(fit, tag) {
  for (cf in colnames(coef(fit))) {
    if (cf == "(Intercept)") next
    t <- topTable(fit, coef=cf, number=Inf, sort.by="none")
    n5  <- sum(t$adj.P.Val < 0.05, na.rm=TRUE)
    n5f <- sum(t$adj.P.Val < 0.05 & abs(t$logFC) > 1, na.rm=TRUE)
    cat(sprintf("  %-28s padj<0.05: %6d   +|LFC|>1: %6d\n", cf, n5, n5f))
    write.csv(t[order(t$adj.P.Val), ],
      paste0(W,"/results/limma_", tag, "_", gsub("[^A-Za-z0-9]","_",cf), ".csv"))
  }
}

des <- model.matrix(~ treatment + egg_mass + month, data=meta)

cat("\n=== B_naive: no blocking, v1-style ===\n")
v <- voom(d, des)
report(eBayes(lmFit(v, des)), "B_naive")

cat("\n=== A_blocked: animal as block, correct for repeated sampling ===\n")
v <- voom(d, des)
cf <- duplicateCorrelation(v, des, block=meta$animal)
cat("consensus within-animal correlation:", round(cf$consensus, 4), "\n")
v <- voom(d, des, block=meta$animal, correlation=cf$consensus)
cf <- duplicateCorrelation(v, des, block=meta$animal)
cat("refined correlation:", round(cf$consensus, 4), "\n")
report(eBayes(lmFit(v, des, block=meta$animal, correlation=cf$consensus)), "A_blocked")

cat("\n=== C_EM3: STP1710.7 vs control within egg mass 3 ===\n")
m3 <- meta[meta$egg_mass=="3" & meta$treatment %in% c("control","STP1710.7"), ]
m3$treatment <- droplevels(m3$treatment); m3$month <- droplevels(m3$month)
d3 <- DGEList(counts=as.matrix(cnt[, m3$stem, drop=FALSE]))
d3 <- calcNormFactors(d3, method="TMM")
des3 <- model.matrix(~ treatment + month, data=m3)
cat("samples:", nrow(m3), " animals:", nlevels(droplevels(m3$animal)), "\n")
v3 <- voom(d3, des3)
c3 <- duplicateCorrelation(v3, des3, block=m3$animal)
v3 <- voom(d3, des3, block=m3$animal, correlation=c3$consensus)
c3 <- duplicateCorrelation(v3, des3, block=m3$animal)
cat("correlation:", round(c3$consensus, 4), "\n")
report(eBayes(lmFit(v3, des3, block=m3$animal, correlation=c3$consensus)), "C_EM3")
cat("\nDONE step17g\n")
