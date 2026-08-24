# Step 17f: one LinDA model per job. Args: tag, model (A|B|C).
# A: treatment + egg_mass + month, (1|animal)   correct for repeated sampling
# B: treatment + egg_mass + month               v1-style, no random effect
# C: STP1710.7 vs control within EM3, (1|animal)
suppressMessages(library(MicrobiomeStat))
a <- commandArgs(trailingOnly=TRUE); tag <- a[1]; mod <- a[2]
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
MINPREV <- 35

meta <- read.delim(paste0(W,"/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
meta$treatment <- relevel(factor(meta$treatment), ref="control")
meta$egg_mass <- factor(meta$egg_mass); meta$month <- factor(meta$month)
meta$animal <- factor(meta$animal)

cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)
cnt <- cnt[, meta$stem, drop=FALSE]
keep <- rowSums(cnt > 0) >= MINPREV
cnt <- cnt[keep, , drop=FALSE]
cat("genes present in >=", MINPREV, "of 44:", nrow(cnt), "\n")

if (mod == "C") {
  meta <- meta[meta$egg_mass=="3" & meta$treatment %in% c("control","STP1710.7"), ]
  meta$treatment <- droplevels(meta$treatment)
  cnt <- cnt[, meta$stem, drop=FALSE]
  f <- "~ treatment + month + (1|animal)"
} else if (mod == "A") {
  f <- "~ treatment + egg_mass + month + (1|animal)"
} else {
  f <- "~ treatment + egg_mass + month"
}
cat("samples:", ncol(cnt), " formula:", f, "\n")

t0 <- Sys.time()
r <- linda(feature.dat=cnt, meta.dat=meta, formula=f,
           feature.dat.type="count", prev.filter=0, mean.abund.filter=0,
           is.winsor=FALSE, n.cores=16, verbose=FALSE)
cat("elapsed:", format(Sys.time()-t0), "\n")
for (v in names(r$output)) {
  o <- r$output[[v]]
  cat(sprintf("  %-28s padj<0.05: %6d   +|LFC|>1: %6d\n", v,
    sum(o$padj < 0.05, na.rm=TRUE),
    sum(o$padj < 0.05 & abs(o$log2FoldChange) > 1, na.rm=TRUE)))
  write.csv(o[order(o$padj), ],
    paste0(W,"/results/linda_", tag, "_", gsub("[^A-Za-z0-9]","_",v), ".csv"))
}
cat("DONE", tag, "\n")
