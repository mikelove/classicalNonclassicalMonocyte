library(DESeq2)
load("txi.rda")
coldata <- y$samples[,"condition",drop=FALSE]
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromTximport(txi, coldata, ~condition)

outlier <- c(12,20) # identified with PCA
keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds <- dds[keep,-outlier]

idx <- c(sample(which(dds$condition == "classical"), 3),
         sample(which(dds$condition == "nonclassical"), 3))
dds$condition[idx]

dds1 <- dds[,idx]
keep2 <- rowSums(counts(dds1) >= 10) >= 3 # to be comparable with edgeR filtering
dds1 <- dds1[keep2,]
dds1 <- DESeq(dds1)
res1 <- results(dds1, name="condition_nonclassical_vs_classical", cooksCutoff=FALSE, independentFiltering=FALSE)

dds2 <- DESeq(dds[keep2,-idx], minReplicates=Inf)
res2 <- results(dds2, name="condition_nonclassical_vs_classical", cooksCutoff=FALSE, independentFiltering=FALSE)

summary(res1)
summary(res2)

table(test=res1$padj < .1, gold=res2$padj < .2)

table(test=res1$padj < .1, gold=sign(res1$log2FoldChange) == sign(res2$log2FoldChange))

# I consider doubling the wrong signs to give FP count, e.g. see this plot:
plot(res1$log2FoldChange[res1$padj < .1], res2$log2FoldChange[res1$padj < .1], cex=.1);
abline(0,1); abline(h=0, v=0)
