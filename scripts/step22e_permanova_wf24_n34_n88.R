# n88 copy of step22e_permanova_wf24_n34.R: 88-metagenome matrix and metadata;
# UHM590.41012 still excluded, so WF24 stays at 34 animals.
# Step 22e: WF24 community tests on 34 animals (exclusions in
# metadata/wf24_excluded.tsv). Abundance as step22a section 4,
# presence/absence as step22c section D.
suppressMessages({library(vegan); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
set.seed(1)
meta <- read.delim(paste0(W,"/metadata/samples_88_treatment.tsv"), stringsAsFactors=FALSE)
ex <- read.delim(paste0(W,"/metadata/wf24_excluded.tsv"), stringsAsFactors=FALSE)
stopifnot(ncol(ex) == 5, nrow(ex) == 2, !("UHM586.41010" %in% meta$stem),
          sum(ex$stem %in% meta$stem) == 1)
m24 <- meta[meta$cohort=="WF24" & !(meta$stem %in% ex$stem), ]
rownames(m24) <- m24$stem
stopifnot(nrow(m24) == 34)
m24$treatment <- factor(m24$treatment)
cat("WF24 groups n =", table(m24$treatment), "\n")
c24 <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary_n88.tsv.gz")),
                  row.names=1, check.names=FALSE)[, m24$stem, drop=FALSE]
cat("\n=== abundance, Bray-Curtis on Hellinger TMM-CPM, genes in >= 30 ===\n")
a24 <- c24[rowSums(c24 > 0) >= 30, , drop=FALSE]
cat("genes:", nrow(a24), " samples:", ncol(a24), "\n")
x24 <- decostand(t(cpm(calcNormFactors(DGEList(counts=as.matrix(a24)), method="TMM"))), "hellinger")
print(adonis2(vegdist(x24, method="bray") ~ treatment, data=m24, permutations=999))
cat("\n=== presence/absence, Jaccard, present at >= 5 reads ===\n")
p24 <- t((c24 >= 5) * 1)
p24 <- p24[, colSums(p24) > 0 & colSums(p24) < nrow(p24), drop=FALSE]
cat("genes variable in presence:", ncol(p24), "\n")
print(adonis2(vegdist(p24, method="jaccard", binary=TRUE) ~ treatment,
              data=m24, permutations=999))
cat("\nDONE step22e_n88\n")
