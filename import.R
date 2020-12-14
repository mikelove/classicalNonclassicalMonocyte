#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF("gencode.v10.annotation.gtf.gz")
#saveDb(txdb, file="gencode.v10.sqlite")
#tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")
#write.table(tx2gene, file="tx2gene.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
tx2gene <- read.table("tx2gene.tsv")

samples <- read.csv("SraRunTable.txt")
samples$condition <- sub("[0-9]+","",samples$Isolate)
samples <- samples[order(samples$condition),]

files <- file.path("quants",samples$Run,"quant.sf")
names(files) <- samples$Run
file.exists(files)

library(tximport)
txi <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")
save(txi, samples, file="txi.rda")

orig.geo.mean <- exp(mean(log(colSums(txi$counts))))

library(edgeR)
y <- DGEList(counts=txi$counts, samples=samples[,c("Run","condition")])
y <- calcNormFactors(y)
edger.cpm <- cpm(y, normalized.lib.sizes=TRUE)

post.geo.mean <- exp(mean(log(colSums(edger.cpm)))) # 1e6 by construction
norm.cts <- round(edger.cpm * orig.geo.mean / post.geo.mean)
write.table(norm.cts, file="norm_cts.tsv", sep="\t", quote=FALSE)

library(DESeq2)
coldata <- y$samples[,"condition",drop=FALSE]
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(norm.cts, colData=coldata, ~condition)
save(dds, file="dds.rda")
