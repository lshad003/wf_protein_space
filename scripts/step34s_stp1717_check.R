# Step 34s check: WF22 STP1717.1 vs control (A_blocked) on c_id95 strict_final, which calls 459 genes (0 on old matrix).
# Model, filter and normalisation as step34e_hostfactor.R WF22 A_blocked: TMM, voom ~ treatment + egg_mass + month,
#   duplicateCorrelation on animal (two rounds). Gene set fixed to the 151,755 tested in the full run (>= 35 of 44)
#   for every refit, so counts are comparable; TMM recomputed on the samples kept.
# Outputs in results/step34s/stp1717_check/: loao_counts.tsv (full refit, leave each STP1717.1 and control animal
#   out, months 2-3 only), m23_STP1717_1.csv.gz (months 2-3 table), refit_full_STP1717_1.csv.gz.
# Optional arguments: fit names to run (full, drop_<animal>, months_2_3); rows then go to loao_counts_<first name>.tsv.
#   Job 30219293 (-p short) timed out after full and 5 drops; the rest were rerun with arguments on epyc.
suppressMessages({library(edgeR); library(limma)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
O <- paste0(W, "/results/step34s/stp1717_check/")
CF <- paste0(W, "/results/step34s/c_id95/step34e_counts_strict_final.tsv.gz")
meta <- read.delim(paste0(W, "/metadata/wf22_design.tsv"), stringsAsFactors=FALSE); rownames(meta) <- meta$stem
cnt <- read.delim(gzfile(CF), row.names=1, check.names=FALSE)[, meta$stem]
a22 <- cnt[rowSums(cnt > 0) >= 35, , drop=FALSE]
cat("genes tested:", nrow(a22), "\n")
COEF <- "treatmentSTP1717.1"
ARGS <- commandArgs(TRUE); DO <- function(n) length(ARGS) == 0 || n %in% ARGS
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
  r <- data.frame(run=tag, samples=nrow(m), animals=nlevels(m$animal), STP1717_animals=length(unique(m$animal[m$treatment=="STP1717.1"])),
    control_animals=length(unique(m$animal[m$treatment=="control"])), correlation=round(cf$consensus, 4),
    padj_lt_0.05=sum(s), up=sum(s & t$logFC > 0), down=sum(s & t$logFC < 0))
  print(r); list(r=r, t=t)
}
R <- list()
if (DO("full")) { x <- fit1(meta, "full"); R[[1]] <- x$r
  write.csv(x$t[order(x$t$adj.P.Val), ], gzfile(paste0(O, "refit_full_STP1717_1.csv.gz"))) }
for (a in sort(unique(meta$animal[meta$treatment %in% c("STP1717.1", "control")]))) {
  g <- unique(meta$treatment[meta$animal == a]); if (!DO(paste0("drop_", a))) next
  R[[length(R) + 1]] <- fit1(meta[meta$animal != a, ], paste0("drop_", a, "_", g))$r
}
if (DO("months_2_3")) { x <- fit1(meta[meta$month %in% c(2, 3), ], "months_2_3"); R[[length(R) + 1]] <- x$r
  write.csv(x$t[order(x$t$adj.P.Val), ], gzfile(paste0(O, "m23_STP1717_1.csv.gz"))) }
OF <- if (length(ARGS) == 0) "loao_counts.tsv" else paste0("loao_counts_", ARGS[1], ".tsv")
write.table(do.call(rbind, R), paste0(O, OF), sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE\n")
