# Step 34s check, array form of scripts/step34s_stp1717_check.R (job 30219293 ran on short and hit the 2 h limit
#   after 4 of 11 refits). Same model: TMM, voom ~ treatment + egg_mass + month, duplicateCorrelation on animal
#   (two rounds), gene set fixed to the 151,755 tested in the full run. One refit per task:
#   0 = full, 1-10 = drop each STP1717.1 / control animal (sorted), 11 = months 2-3 only.
# Output: results/step34s/stp1717_check/loao/task_<i>.tsv (one row); task 11 also writes m23_STP1717_1.csv.gz.
suppressMessages({library(edgeR); library(limma)})
i <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
O <- paste0(W, "/results/step34s/stp1717_check/")
CF <- paste0(W, "/results/step34s/c_id95/step34e_counts_strict_final.tsv.gz")
meta <- read.delim(paste0(W, "/metadata/wf22_design.tsv"), stringsAsFactors=FALSE); rownames(meta) <- meta$stem
cnt <- read.delim(gzfile(CF), row.names=1, check.names=FALSE)[, meta$stem]
a22 <- cnt[rowSums(cnt > 0) >= 35, , drop=FALSE]
cat("genes tested:", nrow(a22), "\n")
COEF <- "treatmentSTP1717.1"
fit1 <- function(m, tag) {
  for (v in c("treatment", "egg_mass", "month", "animal")) m[[v]] <- droplevels(factor(m[[v]]))
  m$treatment <- relevel(m$treatment, ref="control")
  des <- model.matrix(~ treatment + egg_mass + month, data=m)
  if (!is.fullrank(des)) stop("design not full rank: ", tag)
  d <- calcNormFactors(DGEList(counts=as.matrix(a22[, m$stem])), method="TMM")
  v <- voom(d, des); cf <- duplicateCorrelation(v, des, block=m$animal)
  v <- voom(d, des, block=m$animal, correlation=cf$consensus); cf <- duplicateCorrelation(v, des, block=m$animal)
  t <- topTable(eBayes(lmFit(v, des, block=m$animal, correlation=cf$consensus)), coef=COEF, number=Inf, sort.by="none")
  s <- t$adj.P.Val < 0.05
  r <- data.frame(task=i, run=tag, samples=nrow(m), animals=nlevels(m$animal),
    STP1717_animals=length(unique(m$animal[m$treatment=="STP1717.1"])),
    control_animals=length(unique(m$animal[m$treatment=="control"])), correlation=round(cf$consensus, 4),
    padj_lt_0.05=sum(s), up=sum(s & t$logFC > 0), down=sum(s & t$logFC < 0))
  print(r); list(r=r, t=t)
}
an <- sort(unique(meta$animal[meta$treatment %in% c("STP1717.1", "control")]))
if (i == 0) {
  x <- fit1(meta, "full")
} else if (i <= length(an)) {
  a <- an[i]; x <- fit1(meta[meta$animal != a, ], paste0("drop_", a, "_", unique(meta$treatment[meta$animal == a])))
} else {
  x <- fit1(meta[meta$month %in% c(2, 3), ], "months_2_3")
  write.csv(x$t[order(x$t$adj.P.Val), ], gzfile(paste0(O, "m23_STP1717_1.csv.gz")))
}
write.table(x$r, paste0(O, "loao/task_", i, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE\n")
