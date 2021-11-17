library(DESeq2)
load("tximeta_output.rda")
gse$condition <- factor(gse$condition)
dds <- DESeqDataSet(gse, ~condition)

keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds <- dds[keep,]

vsd <- vst(dds)
plotPCA(vsd, ntop=2000)
dds <- dds[,-c(12,20)] # large outliers
vsd <- vst(dds)
plotPCA(vsd, ntop=2000)

# still heterogeneity in PC2, remove unwanted variation
dds_all <- DESeq(dds, minRep=Inf, quiet=TRUE)
res_all <- results(dds_all, cooks=FALSE)
library(RUVSeq)
set <- newSeqExpressionSet(counts(dds))
set <- betweenLaneNormalization(set, which="upper")
not_sig <- rownames(res_all)[which(res_all$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not_sig ]
set <- RUVg(set, empirical, k=2)

for (i in 1:2) {
  w <- paste0("W_",i)
  colData(vsd)[w] <- pData(set)[w]
  colData(dds)[w] <- pData(set)[w]
}
design(dds) <- ~W_1 + W_2 + condition

pca1 <- DESeq2::plotPCA(vsd, intgroup="W_1", ntop=2000)
pca2 <- DESeq2::plotPCA(vsd, intgroup="W_2", ntop=2000)
library(patchwork)
pca1 + pca2

# select 'n' samples from each group
n <- 5
idx <- c(sample(which(dds$condition == "classical"), n),
         sample(which(dds$condition == "nonclassical"), n))
dds$condition[idx]
# small subset samples
dds1 <- dds[,idx]
keep2 <- rowSums(counts(dds1) >= 10) >= n # to be comparable with edgeR filtering
dds1 <- dds1[keep2,]
dds1 <- DESeq(dds1)
res1 <- results(dds1, cooks=FALSE)
res1$padj[is.na(res1$padj)] <- 1
# larger, held-out samples
dds2 <- DESeq(dds[keep2,-idx], minRep=Inf)
res2 <- results(dds2, cooks=FALSE)
summary(res1)
summary(res2)
# test at nominal FDR and see how many sign changes
alpha <- 0.01
sig <- res1$padj < alpha
table(gold=sign(res1$log2FoldChange) == sign(res2$log2FoldChange), test=sig)
plot(res1$log2FoldChange[sig], res2$log2FoldChange[sig], cex=.1)
abline(0,1); abline(h=0, v=0)
