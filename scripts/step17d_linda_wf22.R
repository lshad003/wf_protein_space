# Step 17d: LinDA on WF22, supported genes, three models.
# A: treatment + egg_mass + month, random (1|animal)   <- correct
# B: treatment + egg_mass + month, no random effect    <- v1-style
# C: STP1710.7 vs control within EM3 only, random (1|animal)
suppressMessages(library(MicrobiomeStat))
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"

meta <- read.delim(paste0(W,"/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
meta$treatment <- relevel(factor(meta$treatment), ref="control")
meta$egg_mass  <- factor(meta$egg_mass)
meta$month     <- factor(meta$month)
meta$animal    <- factor(meta$animal)
cat("samples:", nrow(meta), " animals:", nlevels(meta$animal), "\n")

cat("reading count matrix...\n")
cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)
cnt <- cnt[, meta$stem, drop=FALSE]
cat("genes:", nrow(cnt), " cols:", ncol(cnt), "\n")

# prevalence filter: present in at least half the WF22 samples
keep <- rowSums(cnt > 0) >= 22
cnt <- cnt[keep, , drop=FALSE]
cat("genes after prevalence filter (>=22 of 44):", nrow(cnt), "\n")

run <- function(tag, cn, mt, formula) {
  cat("\n=== ", tag, " ===\n", sep="")
  r <- linda(feature.dat=cn, meta.dat=mt, formula=formula,
             feature.dat.type="count", prev.filter=0, mean.abund.filter=0,
             is.winsor=TRUE, outlier.pct=0.03, n.cores=4, verbose=FALSE)
  for (v in names(r$output)) {
    o <- r$output[[v]]
    n5 <- sum(o$padj < 0.05, na.rm=TRUE)
    n5f <- sum(o$padj < 0.05 & abs(o$log2FoldChange) > 1, na.rm=TRUE)
    cat(sprintf("  %-28s padj<0.05: %6d   +|LFC|>1: %6d\n", v, n5, n5f))
    write.csv(o[order(o$padj), ],
      paste0(W,"/results/linda_", tag, "_", gsub("[^A-Za-z0-9]","_",v), ".csv"))
  }
  invisible(r)
}

run("A_mixed", cnt, meta, "~ treatment + egg_mass + month + (1|animal)")
run("B_naive", cnt, meta, "~ treatment + egg_mass + month")

m3 <- meta[meta$egg_mass=="3" & meta$treatment %in% c("control","STP1710.7"), ]
m3$treatment <- droplevels(m3$treatment)
cat("\nEM3 subset: ", nrow(m3), " samples, ", nlevels(droplevels(m3$animal)), " animals\n", sep="")
run("C_EM3", cnt[, m3$stem, drop=FALSE], m3, "~ treatment + month + (1|animal)")
cat("\nDONE step17d\n")
