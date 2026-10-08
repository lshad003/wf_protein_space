# Step 17i n88: copy of step17i_limma_wf24_n34.R on the 88-metagenome matrix
# and supported set. UHM586.41010 is already absent from the n88 metadata and
# matrix; UHM590.41012 is excluded through metadata/wf24_excluded.tsv, so WF24
# stays at 34 animals. Same model and prevalence filter as step17i.
suppressMessages({library(limma); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
meta <- read.delim(paste0(W,"/metadata/samples_88_treatment.tsv"), stringsAsFactors=FALSE)
ex <- read.delim(paste0(W,"/metadata/wf24_excluded.tsv"), stringsAsFactors=FALSE)
stopifnot(ncol(ex) == 5, nrow(ex) == 2, !("UHM586.41010" %in% meta$stem),
          sum(ex$stem %in% meta$stem) == 1)
m24 <- meta[meta$cohort=="WF24" & !(meta$stem %in% ex$stem), ]
rownames(m24) <- m24$stem
stopifnot(nrow(m24) == 34)
m24$treatment <- relevel(factor(m24$treatment), ref="4")
cat("excluded:", ex$stem, "\n")
cat("WF24 groups:", levels(m24$treatment), "n =", table(m24$treatment), "\n")
cat("WF24 egg masses n =", table(m24$egg_mass), "\n")
cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary_n88.tsv.gz")),
                  row.names=1, check.names=FALSE)[, m24$stem, drop=FALSE]
cnt <- cnt[rowSums(cnt > 0) >= 30, , drop=FALSE]
cat("samples:", ncol(cnt), " genes present in >= 30:", nrow(cnt), "\n")
d <- calcNormFactors(DGEList(counts=as.matrix(cnt)), method="TMM")
des <- model.matrix(~ treatment, data=m24)
fit <- eBayes(lmFit(voom(d, des), des))
for (cf in colnames(coef(fit))[-1]) {
  t <- topTable(fit, coef=cf, number=Inf, sort.by="none")
  cat(sprintf("  %-14s padj<0.05: %6d   +|LFC|>1: %6d\n", cf,
    sum(t$adj.P.Val < 0.05, na.rm=TRUE),
    sum(t$adj.P.Val < 0.05 & abs(t$logFC) > 1, na.rm=TRUE)))
  write.csv(t[order(t$adj.P.Val), ],
    paste0(W,"/results/limma_WF24_n34_", gsub("[^A-Za-z0-9]","_",cf), "_n88.csv"))
}
cat("\nDONE step17i_n88\n")
